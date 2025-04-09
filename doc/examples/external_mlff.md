# Advanced Usage: Training MLFF on-the-fly during metadynamics simulation

The socket makes porting SPARC to external codes with python interfaces trivial. The [PLUMED](https://github.com/plumed/plumed2) package's existing `ASE` interface can be coupled with the SPARC-X-API to train MLFF on-the-fly during metadynamics simulations.

The ASE interface to PLUMED requires both the PLUMED binary and
`py-plumed` library to work properly, check [ASE-PLUMED
documentation](https://wiki.fysik.dtu.dk/ase//ase/calculators/plumed.html)
for more details. The most convenient way to install the extra
dependencies is through conda:
```{code} bash
conda install -c conda-forge py-plumed
```

```{note}
Check the [PLUMED documentation](https://www.plumed.org/doc) for other installation options. In addition, the environmental variable `PLUMED_KERNEL` must be set to the location of the shared library `libplumedKernel.so` to work properly.
```


The setup on the SPARC side remains unchanged from other socket examples:

```python
from sparc.calculator import SPARC
from ase import units
from ase import Atoms
from ase.calculators.plumed import Plumed
from ase.md.nvtberendsen import NVTBerendsen

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1,1,1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "MLFF_FLAG": 1,
    "MLFF_INITIAL_STEPS_TRAIN": 5,
}

# MD parameters
timestep = 2 * units.fs  # the units module contains the conversion factor to go from units.<unit> to ASE units via multiplication
ps = 1000 * units.fs
```

Here, we will consider a 5 atom Ag cluster.

```python
Ag_cluster = Atoms('Ag5', positions = [(0.0, 2.6579, 0.9366), (0.0, -1.3587, -1.4045), (0.0, 0.0, 0.9358), (0.0, -2.6579, 0.9366),
                                      (0.0, 1.3587, -1.4045)],
                  pbc = (0,0,0))
Ag_cluster.set_cell([20., 24., 24.])
Ag_cluster.center()
```

We include PLUMED specific parameters to initialize the metadynamics simulations, taking care to be consistent with native `ASE` and PLUMED units.


```python
setup = [
    f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",  # Set units to match desired properties

    # Calculate the center of mass of atoms 1-5
    "com: COM ATOMS=1-5",

    # Define the coordination number (C)
    "c: COORDINATION GROUPA=1-5 SWITCH={RATIONAL R_0=3.0 NN=8 MM=16} NOPBC",

    # Define the radius of gyration (R)
    "r: GYRATION TYPE=RADIUS ATOMS=1-5 NOPBC",

    # Compute CV1 and CV2 as orthogonal linear combinations of C and R
    "cv1: COMBINE ARG=c,r COEFFICIENTS=0.99715,-0.07534 PERIODIC=NO",
    "cv2: COMBINE ARG=c,r COEFFICIENTS=0.07534,0.99715 PERIODIC=NO",

    # Apply lower wall on CV1 at 5.0 with harmonic constant 10 eV
    f"LOWER_WALLS ARG=cv1 AT=5.0 KAPPA=10.0",

    # Apply lower wall on CV2 at 3.0 with harmonic constant 50 eV
    f"UPPER_WALLS ARG=cv2 AT=3.0 KAPPA=50.0",

    # Perform well-tempered metadynamics on CV1 and CV2
    f"METAD ARG=cv1,cv2 HEIGHT=0.3 PACE=500 SIGMA=0.3,0.03 GRID_MIN=0.0,0.0 GRID_MAX=10.0,5.0 GRID_BIN=500,500 BIASFACTOR=100 FILE=HILLS",

    # Print out the collective variables for monitoring
    "PRINT ARG=cv1,cv2,c,r FILE=COLVAR STRIDE=1"
]
```

And finally we run a brief simulation in socket mode:

```python
with SPARC(use_socket=True, **calc_params) as calc:
    Ag_cluster.calc = Plumed(calc=calc,
                    input=setup,
                    timestep=timestep,
                    atoms=Ag_cluster,
                    kT=0.00861733) # 10 K in eV thermal energy units
    dyn = NVTBerendsen(Ag_cluster, timestep, temperature_K=0.00861733/units.kB, taut=50*units.fs,
                fixcm=False, trajectory='Ag-cluster-metadynamics.traj')
    dyn.run(10)
```

The full script can be downloaded [here](https://raw.githubusercontent.com/SPARC-X/SPARC-X-API/refs/heads/master/examples/socket/md/mlff/plumed/run.py).
