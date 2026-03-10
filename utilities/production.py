"""
production.py
=============
Production MD stage for membrane protein simulations using OpenMM.
Runs N sequential replicas, each in its own subdirectory.

Expected directory layout (produced by equilibration.py)
---------------------------------------------------------
::

    run_dir/
    ├── build/
    │   └── structure.psf            ← topology
    ├── common/toppar/               ← CHARMM36 parameter files
    └── equilibration/
        ├── equil_boxsize.csv        ← box dimensions
        └── equil_coord.pickle       ← final equilibrated coordinates

After a successful run, outputs are written to::

    run_dir/
    └── production/
        ├── rep1/
        │   ├── prod_traj.dcd        ← trajectory (unwrapped)
        │   ├── production.log       ← thermodynamic state data
        │   ├── checkpoint_vel.pickle
        │   ├── checkpoint_coor.pickle
        │   └── last_frame.pdb       ← quick visual check
        ├── rep2/
        │   └── …
        └── repN/
            └── …

Usage
-----
::

    python production.py <run_dir> <start|continue> [--replicas N]

Arguments
---------
run_dir        : str
    Root simulation directory (see layout above).
mode           : "start" | "continue"
    ``start``    – draw velocities from Maxwell-Boltzmann; begin fresh DCD.
    ``continue`` – reload positions/velocities from checkpoints; append to DCD.
--replicas N   : int  (optional, default = 5)
    Number of sequential replicas to run (rep1 … repN).

Notes
-----
* Hydrogen Mass Repartitioning (HMR, 4 amu) + AllBonds constraints → 4 fs step.
* Replicas share the same starting coordinates (equil_coord.pickle) but
  receive independently drawn velocities, producing statistically independent
  trajectories.
* In ``continue`` mode every replica resumes from its own checkpoint files.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import os
import sys
import csv
import pickle
import argparse

from openmm.app import (
    CharmmPsfFile, CharmmParameterSet, PDBFile,
    PME, AllBonds, Simulation, StateDataReporter, DCDReporter,
)
from openmm import LangevinIntegrator
import openmm.unit as u

# ---------------------------------------------------------------------------
# Simulation parameters  –  edit here if needed
# ---------------------------------------------------------------------------
PROD_STEPS       = 12_500_000   # steps per replica  (~50 ns at 4 fs)
SAVE_EVERY       =     50_000   # steps between DCD frames
LOG_EVERY        =      5_000   # steps between log entries
TEMPERATURE      = 310 * u.kelvin
TIMESTEP         = 4 * u.femtoseconds
FRICTION         = 0.1 / u.picoseconds
NONBONDED_CUTOFF = 9.0 * u.angstroms
SWITCH_DISTANCE  = 7.5 * u.angstroms
HYDROGEN_MASS    = 4 * u.amu
N_REPLICAS       = 5   # default number of sequential replicas

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Run sequential production replicas for a membrane protein simulation.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=(
        "Examples:\n"
        "  python production.py /path/to/run_dir start\n"
        "  python production.py /path/to/run_dir start --replicas 3\n"
        "  python production.py /path/to/run_dir continue --replicas 5\n"
    ),
)
parser.add_argument("run_dir", help="Root simulation directory.")
parser.add_argument(
    "mode",
    choices=["start", "continue"],
    help="'start' = fresh run; 'continue' = resume from checkpoints.",
)
parser.add_argument(
    "--replicas",
    type=int,
    default=N_REPLICAS,
    metavar="N",
    help=f"Number of sequential replicas to run (default: {N_REPLICAS}).",
)

args     = parser.parse_args()
run_dir  = os.path.realpath(args.run_dir)
mode     = args.mode
n_reps   = args.replicas
is_start = mode == "start"

if n_reps < 1:
    sys.exit("ERROR: --replicas must be a positive integer.")

# ---------------------------------------------------------------------------
# Derived shared paths
# ---------------------------------------------------------------------------
build_dir  = os.path.join(run_dir, "build")
equil_dir  = os.path.join(run_dir, "equilibration")
toppar_dir = os.path.join(run_dir, "common", "toppar")
prod_dir   = os.path.join(run_dir, "production")

psf_path      = os.path.join(build_dir,  "structure.psf")
equil_boxsize = os.path.join(equil_dir,  "equil_boxsize.csv")
equil_coords  = os.path.join(equil_dir,  "equil_coord.pickle")

toppar_files = [
    os.path.join(toppar_dir, "toppar_water_ions.str"),
    os.path.join(toppar_dir, "top_all36_prot.rtf"),
    os.path.join(toppar_dir, "top_all36_lipid.rtf"),
    os.path.join(toppar_dir, "par_all36_lipid.prm"),
    os.path.join(toppar_dir, "par_all36_prot.prm"),
]

# ---------------------------------------------------------------------------
# Input validation  –  fail fast before loading OpenMM
# ---------------------------------------------------------------------------
errors = []

if not os.path.isdir(run_dir):
    errors.append(
        f"ERROR: run_dir does not exist or is not a directory:\n"
        f"         {run_dir}\n"
        f"       Run equilibration.py first to prepare this directory."
    )

if os.path.isdir(run_dir):
    # Files always required
    required = {
        psf_path      : "topology from the building pipeline",
        equil_boxsize : "box dimensions written by equilibration.py",
        equil_coords  : "final coordinates written by equilibration.py",
    }
    for fpath, description in required.items():
        if not os.path.isfile(fpath):
            errors.append(
                f"ERROR: required file not found: {fpath}\n"
                f"       Expected: {description}"
            )

    # In continue mode, every replica must have its checkpoint files
    if not is_start:
        for rep in range(1, n_reps + 1):
            rep_dir = os.path.join(prod_dir, f"rep{rep}")
            for fname in ("checkpoint_coor.pickle", "checkpoint_vel.pickle"):
                fpath = os.path.join(rep_dir, fname)
                if not os.path.isfile(fpath):
                    errors.append(
                        f"ERROR: checkpoint not found for rep{rep}:\n"
                        f"         {fpath}\n"
                        f"       Run in 'start' mode first to create checkpoints."
                    )

if errors:
    print("\n".join(errors), file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def read_boxsize(csv_path: str) -> list:
    """
    Read periodic box dimensions from a CSV file (header: x,y,z).

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


def read_pickle(path: str):
    """Load and return a pickled object (positions or velocities)."""
    with open(path, "rb") as fh:
        return pickle.load(fh)


def write_pickle(obj, path: str) -> None:
    """Pickle ``obj`` to ``path``."""
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)


def run_replica(rep_id: int, psf, params, is_start: bool) -> None:
    """
    Build and run a single production replica.

    The replica reads shared inputs (box size, equilibrated coordinates or
    previous checkpoints) and writes all outputs to its own subdirectory
    ``production/rep<rep_id>/``.

    Parameters
    ----------
    rep_id   : int   – replica number (1-based)
    psf      : openmm.app.CharmmPsfFile  – pre-loaded topology
    params   : openmm.app.CharmmParameterSet
    is_start : bool  – True = fresh start; False = continue from checkpoint
    """
    rep_dir = os.path.join(prod_dir, f"rep{rep_id}")
    os.makedirs(rep_dir, exist_ok=True)

    # Per-replica file paths
    out_dcd    = os.path.join(rep_dir, "prod_traj.dcd")
    out_log    = os.path.join(rep_dir, "production.log")
    out_pdb    = os.path.join(rep_dir, "last_frame.pdb")
    ckpt_vel   = os.path.join(rep_dir, "checkpoint_vel.pickle")
    ckpt_coor  = os.path.join(rep_dir, "checkpoint_coor.pickle")

    print(f"\n{'─' * 60}")
    print(f"  Replica {rep_id}/{n_reps}   ({'start' if is_start else 'continue'})")
    print(f"  Output dir : {rep_dir}")
    print(f"{'─' * 60}")

    # --- Set box from equilibration output ---
    psf.setBox(*read_boxsize(equil_boxsize))

    # --- Build OpenMM System with HMR (enables 4 fs time-step) ---
    system = psf.createSystem(
        params,
        nonbondedMethod = PME,
        nonbondedCutoff = NONBONDED_CUTOFF,
        switchDistance  = SWITCH_DISTANCE,
        constraints     = AllBonds,       # required with HMR
        hydrogenMass    = HYDROGEN_MASS,  # hydrogen mass repartitioning
        rigidWater      = True,
    )

    integrator = LangevinIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(psf.topology, system, integrator)

    # --- Set initial positions and velocities ---
    if is_start:
        print("  Loading equilibrated coordinates and drawing velocities …")
        simulation.context.setPositions(read_pickle(equil_coords))
        simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    else:
        print("  Loading checkpoint positions and velocities …")
        simulation.context.setPositions(read_pickle(ckpt_coor))
        simulation.context.setVelocities(read_pickle(ckpt_vel))

    # --- Attach reporters ---
    simulation.reporters.append(
        DCDReporter(
            out_dcd,
            SAVE_EVERY,
            enforcePeriodicBox = False,    # keep molecules whole across frames
            append             = not is_start,
        )
    )
    simulation.reporters.append(
        StateDataReporter(
            out_log, LOG_EVERY,
            step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            progress=True, remainingTime=True, speed=True,
            elapsedTime=True, separator=";", totalSteps=PROD_STEPS,
        )
    )

    # --- Run ---
    print(f"  Running {PROD_STEPS:,} steps …")
    simulation.step(PROD_STEPS)

    # --- Save checkpoints and final PDB ---
    print("  Saving checkpoints and last frame …")
    write_pickle(
        simulation.context.getState(getVelocities=True).getVelocities(),
        ckpt_vel,
    )
    write_pickle(
        simulation.context.getState(getPositions=True).getPositions(),
        ckpt_coor,
    )
    PDBFile.writeFile(
        simulation.topology,
        simulation.context.getState(getPositions=True).getPositions(),
        open(out_pdb, "w"),
    )

    print(f"  Replica {rep_id} complete.")
    print(f"    Trajectory : {out_dcd}")
    print(f"    Log        : {out_log}")
    print(f"    Last frame : {out_pdb}")


# ---------------------------------------------------------------------------
# Main  –  load shared resources once, then run replicas sequentially
# ---------------------------------------------------------------------------

print("=" * 60)
print("  Production")
print("=" * 60)
print(f"  Run directory : {run_dir}")
print(f"  Mode          : {mode}")
print(f"  Replicas      : {n_reps}  (rep1 … rep{n_reps}, sequential)")
print(f"  Steps / rep   : {PROD_STEPS:,}  (~{PROD_STEPS * 4e-6:.0f} ns)")
print()

# Load PSF and force-field parameters once — reused for every replica
# Note: the equilibration PDB is NOT loaded here; OpenMM ≥ 8.0 can produce
# duplicate-atom errors from saved PDB files; coordinates come via pickle.
print("Loading PSF and CHARMM36 parameters …")
psf    = CharmmPsfFile(psf_path)
params = CharmmParameterSet(*toppar_files)

os.makedirs(prod_dir, exist_ok=True)

for rep in range(1, n_reps + 1):
    run_replica(rep, psf, params, is_start)

print()
print("=" * 60)
print(f"  All {n_reps} replica(s) complete.")
print(f"  Results in : {prod_dir}")
print("=" * 60)
