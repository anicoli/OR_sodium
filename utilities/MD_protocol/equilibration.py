"""
equilibration.py
================
Equilibration stage for membrane protein MD simulations using OpenMM.

Expected directory layout
--------------------------
::

    run_dir/
    ├── build/
    │   ├── structure.psf        ← input topology
    │   └── structure.pdb        ← input coordinates
    └── common/
        └── toppar/              ← CHARMM36 parameter files
            ├── toppar_water_ions.str
            ├── top_all36_prot.rtf
            ├── top_all36_lipid.rtf
            ├── par_all36_lipid.prm
            └── par_all36_prot.prm

After a successful run, outputs are written to::

    run_dir/
    └── equilibration/
        ├── equil.xtc              ← trajectory
        ├── equil_boxsize.csv      ← final box dimensions (read by production.py)
        ├── equil_coord.pickle     ← final coordinates   (read by production.py)
        └── equilibration.log      ← thermodynamic state data

Usage
-----
::

    python equilibration.py <run_dir> <boxsize_csv>

Arguments
---------
run_dir      : str
    Root simulation directory (see layout above).
boxsize_csv  : str
    CSV file supplying the initial periodic box dimensions in Angstroms.
    Must contain a header row followed by exactly one data row::

        x,y,z
        90.0,90.0,100.0

Restraint schedule
------------------
Harmonic positional restraints (backbone + internal waters, chain 1) are
applied and then gradually released over the course of the run:

==========  =====================================================
0 – 50 %    Full restraint  (k = RESTRAINT_FORCE)
50 – 75 %   Linear ramp     (k → 0)
75 – 100 %  No restraint    (k = 0)
==========  =====================================================

References
----------
[1] MonteCarloMembraneBarostat surface tension:
    http://docs.openmm.org/latest/userguide/application/02_running_sims.html
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import os
import sys
import csv
import pickle

import mdtraj as md
from openmm.app import (
    CharmmPsfFile, PDBFile, CharmmParameterSet,
    PME, HBonds, Simulation, StateDataReporter,
)
from openmm import (
    System, CustomExternalForce, MonteCarloMembraneBarostat,
    LangevinIntegrator,
)
import openmm.unit as u

# ---------------------------------------------------------------------------
# Simulation parameters  –  edit here if needed
# ---------------------------------------------------------------------------
EQUIL_STEPS      = 20_000_000   # total integration steps
SAVE_EVERY       =     25_000   # steps between trajectory frames
LOG_EVERY        =      5_000   # steps between log entries
TEMPERATURE      = 310 * u.kelvin
PRESSURE         = 1.01325 * u.bar
TIMESTEP         = 2 * u.femtoseconds
FRICTION         = 0.1 / u.picoseconds
NONBONDED_CUTOFF = 9.0 * u.angstroms
SWITCH_DISTANCE  = 7.5 * u.angstroms
RESTRAINT_FORCE  = 1.0   # kcal mol⁻¹ Å⁻²  (initial; released gradually)
RESTRAINT_SEL    = "backbone or (water and chainid 1)"  # chainid 1 = internal waters

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
if len(sys.argv) != 3:
    sys.exit(
        "Usage: python equilibration.py <run_dir> <boxsize_csv>\n"
        "\n"
        "  run_dir     – root simulation directory; must contain:\n"
        "                  build/structure.psf\n"
        "                  build/structure.pdb\n"
        "                  common/toppar/  (CHARMM36 parameter files)\n"
        "  boxsize_csv – CSV with initial box dimensions (Å), e.g.:\n"
        "                  x,y,z\n"
        "                  90.0,90.0,100.0\n"
    )

run_dir     = os.path.realpath(sys.argv[1])
boxsize_csv = os.path.realpath(sys.argv[2])

# ---------------------------------------------------------------------------
# Input validation  –  fail fast with clear messages before loading OpenMM
# ---------------------------------------------------------------------------
errors = []

if not os.path.isdir(run_dir):
    errors.append(
        f"ERROR: run_dir does not exist or is not a directory:\n"
        f"         {run_dir}\n"
        f"       Create it and ensure it contains build/ and common/toppar/."
    )

if not os.path.isfile(boxsize_csv):
    errors.append(
        f"ERROR: boxsize CSV not found:\n"
        f"         {boxsize_csv}\n"
        f"       The file must contain a header 'x,y,z' and one data row, e.g.:\n"
        f"         x,y,z\n"
        f"         90.0,90.0,100.0"
    )

if os.path.isdir(run_dir):
    for fname in ("build/structure.psf", "build/structure.pdb"):
        fpath = os.path.join(run_dir, fname)
        if not os.path.isfile(fpath):
            errors.append(
                f"ERROR: required input file not found:\n"
                f"         {fpath}\n"
                f"       Make sure the building pipeline has completed successfully."
            )

if errors:
    print("\n".join(errors), file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------

# --- Inputs ---
build_dir  = os.path.join(run_dir, "build")
psf_path   = os.path.join(build_dir, "structure.psf")
pdb_path   = os.path.join(build_dir, "structure.pdb")
toppar_dir = os.path.join(run_dir, "common", "toppar")

toppar_files = [
    os.path.join(toppar_dir, "toppar_water_ions.str"),
    os.path.join(toppar_dir, "top_all36_prot.rtf"),
    os.path.join(toppar_dir, "top_all36_lipid.rtf"),
    os.path.join(toppar_dir, "par_all36_lipid.prm"),
    os.path.join(toppar_dir, "par_all36_prot.prm"),
]

# --- Outputs  (all inside run_dir/equilibration/) ---
equil_dir = os.path.join(run_dir, "equilibration")
os.makedirs(equil_dir, exist_ok=True)

out_xtc     = os.path.join(equil_dir, "equil.xtc")
out_boxsize = os.path.join(equil_dir, "equil_boxsize.csv")
out_coords  = os.path.join(equil_dir, "equil_coord.pickle")
out_log     = os.path.join(equil_dir, "equilibration.log")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def read_boxsize(csv_path: str) -> list:
    """
    Read periodic box dimensions from a CSV file.

    The file must have a header row ``x,y,z`` followed by one data row
    with the three box lengths in Angstroms.

    Parameters
    ----------
    csv_path : str

    Returns
    -------
    list of openmm.unit.Quantity  –  [x, y, z] in Angstroms
    """
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        row = next(reader)
        return [
            u.Quantity(float(row["x"]), u.angstrom),
            u.Quantity(float(row["y"]), u.angstrom),
            u.Quantity(float(row["z"]), u.angstrom),
        ]


def add_positional_restraints(
    system: System,
    pdb: PDBFile,
    pdb_path: str,
    selection: str,
    force_constant: float,
) -> None:
    """
    Add harmonic positional restraints to a selected set of atoms.

    Uses a ``CustomExternalForce`` with a periodic-distance potential so
    restraints behave correctly across PBC boundaries.  The global
    parameter ``k`` can be updated at runtime to gradually release them.

    Parameters
    ----------
    system        : openmm.System
    pdb           : openmm.app.PDBFile
    pdb_path      : str    – used by MDTraj to resolve the atom selection
    selection     : str    – MDTraj atom selection string
    force_constant: float  – initial k in kcal mol⁻¹ Å⁻²
    """
    restraint = CustomExternalForce(
        "k * periodicdistance(x, y, z, x0, y0, z0)^2"
    )
    restraint.addGlobalParameter(
        "k", force_constant * u.kilocalories_per_mole / (u.angstrom ** 2)
    )
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    restrained_indices = set(md.load(pdb_path).topology.select(selection))
    for atom in pdb.topology.atoms():
        if atom.index in restrained_indices:
            restraint.addParticle(atom.index, pdb.positions[atom.index])

    system.addForce(restraint)


def write_boxsize(simulation: Simulation, csv_path: str) -> None:
    """
    Write the current periodic box dimensions (Å) to a CSV file.

    Parameters
    ----------
    simulation : openmm.app.Simulation
    csv_path   : str
    """
    box = simulation.context.getState(getParameters=True).getPeriodicBoxVectors()
    x = box[0][0].value_in_unit(u.angstrom)
    y = box[1][1].value_in_unit(u.angstrom)
    z = box[2][2].value_in_unit(u.angstrom)

    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["x", "y", "z"])
        writer.writerow([x, y, z])


def write_coordinates(simulation: Simulation, pickle_path: str) -> None:
    """
    Pickle the current atom positions for use in production.py.

    Parameters
    ----------
    simulation  : openmm.app.Simulation
    pickle_path : str
    """
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(pickle_path, "wb") as fh:
        pickle.dump(positions, fh)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 60)
print("  Equilibration")
print("=" * 60)
print(f"  Run directory  : {run_dir}")
print(f"  Box-size CSV   : {boxsize_csv}")
print(f"  Output dir     : {equil_dir}")
print()

# --- Load structure and force field ---
print("Loading PSF, PDB, and CHARMM36 parameters …")
psf    = CharmmPsfFile(psf_path)
pdb    = PDBFile(pdb_path)
params = CharmmParameterSet(*toppar_files)

# --- Set periodic box from the supplied CSV ---
psf.setBox(*read_boxsize(boxsize_csv))

# --- Build OpenMM System ---
print("Building OpenMM System …")
system = psf.createSystem(
    params,
    nonbondedMethod = PME,
    nonbondedCutoff = NONBONDED_CUTOFF,
    switchDistance  = SWITCH_DISTANCE,
    constraints     = HBonds,
    rigidWater      = True,
)

add_positional_restraints(system, pdb, pdb_path, RESTRAINT_SEL, RESTRAINT_FORCE)

# Monte Carlo membrane barostat: XY isotropic, Z free (membrane normal)  [1]
system.addForce(
    MonteCarloMembraneBarostat(
        PRESSURE,
        0,              # surface tension (dyn/cm); 0 = tension-free
        TEMPERATURE,
        MonteCarloMembraneBarostat.XYIsotropic,
        MonteCarloMembraneBarostat.ZFree,
    )
)

integrator = LangevinIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

# --- Initialise, minimise, attach reporters ---
print("Building Simulation and minimising energy …")
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy(maxIterations=5_000)

simulation.reporters.append(md.reporters.XTCReporter(out_xtc, SAVE_EVERY))
simulation.reporters.append(
    StateDataReporter(
        out_log, LOG_EVERY,
        step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        elapsedTime=True, separator=";", totalSteps=EQUIL_STEPS,
    )
)

# --- Run with gradual restraint release ---
# The run is divided into 100 equal chunks; k is updated each chunk:
#   chunks  0 – 49  : k = RESTRAINT_FORCE  (full restraint)
#   chunks 50 – 75  : k ramps linearly RESTRAINT_FORCE → 0
#   chunks 76 – 99  : k = 0                (fully unrestrained)
print(f"Running equilibration ({EQUIL_STEPS:,} steps) …")
steps_per_chunk = EQUIL_STEPS // 100

for chunk in range(100):
    if chunk < 50:
        k = RESTRAINT_FORCE
    elif chunk <= 75:
        k = ((75 - chunk) / 25.0) * RESTRAINT_FORCE
    else:
        k = 0.0
    simulation.context.setParameter("k", k)
    simulation.step(steps_per_chunk)

# --- Save outputs ---
print("Saving outputs …")
write_boxsize(simulation, out_boxsize)
write_coordinates(simulation, out_coords)

print()
print("Equilibration complete.")
print(f"  Trajectory    : {out_xtc}")
print(f"  Box size      : {out_boxsize}")
print(f"  Coordinates   : {out_coords}")
print(f"  Log           : {out_log}")
