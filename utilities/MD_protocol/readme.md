# MD Protocol — Odorant receptors

Membrane protein molecular dynamics simulations using [OpenMM](https://openmm.org/).  
This folder contains the equilibration and production scripts used to generate the OR–sodium trajectories described in the associated publication.

---

## Dependencies

| Package     | Version tested | Notes                              |
|-------------|----------------|------------------------------------|
| OpenMM      | ≥ 8.0          | GPU-accelerated MD engine          |
| MDTraj      | ≥ 1.9          | Atom selection + XTC reporter      |
| Python      | ≥ 3.9          |                                    |

Install via conda (recommended):

```bash
conda install -c conda-forge openmm mdtraj
```

---

## Input files

Before running, the simulation directory (`run_dir`) must have the following layout, typically produced by the building pipeline (HTMD/VMD/CHARMM-GUI):

```
run_dir/
├── build/
│   ├── structure.psf        # CHARMM topology
│   └── structure.pdb        # Initial coordinates
└── common/
    └── toppar/              # CHARMM36 force-field parameter files
        ├── toppar_water_ions.str
        ├── top_all36_prot.rtf
        ├── top_all36_lipid.rtf
        ├── par_all36_lipid.prm
        └── par_all36_prot.prm
```

You also need a `boxsize.csv` file with the initial periodic box dimensions in Ångströms:

```
x,y,z
90.0,90.0,100.0
```

---

## Workflow

```
build/structure.psf + structure.pdb
            │
            ▼
    equilibration.py          →   run_dir/equilibration/
            │
            ▼
    production.py             →   run_dir/production/rep1/ … repN/
```

### 1 — Equilibration

```bash
python equilibration.py <run_dir> <boxsize_csv>
```

**What it does:**

- Loads `build/structure.psf` and `build/structure.pdb`
- Energy-minimises the system (up to 5 000 steps)
- Runs **40 ns** of NPT equilibration (2 fs time-step, 310 K, 1 atm)  
  with a Monte Carlo membrane barostat (XY isotropic, Z free)
- Applies harmonic positional restraints on backbone atoms and internal  
  waters (chain 1), then releases them gradually:

  | Simulation progress | Restraint k            |
  |---------------------|------------------------|
  | 0 – 50 %            | Full (1 kcal/mol/Å²)   |
  | 50 – 75 %           | Linear ramp → 0        |
  | 75 – 100 %          | None                   |

**Outputs** (written to `run_dir/equilibration/`):

| File                  | Description                                      |
|-----------------------|--------------------------------------------------|
| `equil.xtc`           | Trajectory                                       |
| `equil_boxsize.csv`   | Final box dimensions (passed to production.py)   |
| `equil_coord.pickle`  | Final atom positions (passed to production.py)   |
| `equilibration.log`   | Thermodynamic state data (`;`-separated CSV)     |

---

### 2 — Production

```bash
python production.py <run_dir> start [--replicas N]
```

To resume from checkpoints:

```bash
python production.py <run_dir> continue [--replicas N]
```

| Argument      | Default | Description                                       |
|---------------|---------|---------------------------------------------------|
| `run_dir`     | —       | Root simulation directory                         |
| `start`       | —       | Fresh run; draws velocities from Maxwell–Boltzmann|
| `continue`    | —       | Resumes from per-replica checkpoint files         |
| `--replicas N`| 5       | Number of sequential independent replicas         |

**What it does:**

- Runs **N replicas sequentially** (rep1 → rep2 → … → repN)
- Each replica starts from the same equilibrated coordinates but with  
  independently drawn velocities, producing statistically independent trajectories
- Uses Hydrogen Mass Repartitioning (HMR, 4 amu) with a 4 fs time-step
- Each replica runs for **50 ns** (12 500 000 steps)

**Outputs** (written to `run_dir/production/repN/`):

| File                      | Description                               |
|---------------------------|-------------------------------------------|
| `prod_traj.dcd`           | Trajectory (unwrapped, DCD format)        |
| `production.log`          | Thermodynamic state data                  |
| `last_frame.pdb`          | Final frame (quick visual check)          |
| `checkpoint_coor.pickle`  | Positions checkpoint (for continuation)   |
| `checkpoint_vel.pickle`   | Velocities checkpoint (for continuation)  |

---

## Output directory structure

After a complete run with 5 replicas:

```
run_dir/
├── build/
├── common/
├── equilibration/
│   ├── equil.xtc
│   ├── equil_boxsize.csv
│   ├── equil_coord.pickle
│   └── equilibration.log
└── production/
    ├── rep1/
    │   ├── prod_traj.dcd
    │   ├── production.log
    │   ├── last_frame.pdb
    │   ├── checkpoint_coor.pickle
    │   └── checkpoint_vel.pickle
    ├── rep2/
    ├── rep3/
    ├── rep4/
    └── rep5/
```

---

## Simulation parameters

All key parameters are defined at the top of each script and can be  
edited without touching the simulation logic:

| Parameter         | Equilibration   | Production       |
|-------------------|-----------------|------------------|
| Temperature       | 310 K           | 310 K            |
| Time-step         | 2 fs            | 4 fs (HMR)       |
| Nonbonded cutoff  | 9 Å             | 9 Å              |
| Switch distance   | 7.5 Å           | 7.5 Å            |
| Constraints       | HBonds          | AllBonds         |
| Barostat          | MC membrane NPT | None             |
| Total steps       | 20 000 000      | 12 500 000 / rep |
| Trajectory output | XTC, every 25k  | DCD, every 50k   |

---

## Force field

CHARMM36 force field for proteins and lipids, with TIP3P water and ions.  
Parameter files must be provided by the user (not included in this repository  
for licensing reasons). They can be obtained from the  
[CHARMM-GUI](https://www.charmm-gui.org/) or the  
[MacKerell lab website](https://mackerell.umaryland.edu/charmm_ff.shtml).

> P. Eastman et al. *OpenMM 8: Molecular Dynamics Simulation with Machine  
> Learning Potentials.* J. Phys. Chem. B, 2024.
