# Examples

For template scripts utilizing the API in socket mode and via FileIO mode, see the `examples/` directory in the [github repo](https://github.com/SPARC-X/SPARC-X-API/tree/master/examples).


## Simple DFT workflows with SPARC-X-API

As documentation for [basic usage](basic_usage.md) shows, replacing an
existing DFT-based workflow using SPARC-X-API calculator interface can
be as easy as simply swapping the calculator instance to `sparc.SPARC`
from other codes like VASP, QE or GPAW.

Here we show a simple example that calculates the equation of state
(EOS) of bulk aluminum and determine its optimal lattice constant,
adapted from GPAW's
[tutorial](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/lattice_constants/lattice_constants.html). In
the GPAW tutorial, a GPAW calculator in Planewave (PW) mode is created
like follows:

```{code} python
from gpaw import GPAW, PW
calc = GPAW(mode=PW(ecut),
            xc='PBE',
			kpts=(8, 8, 8),
			basis='dzp',
			txt=f'Al-{ecut}.txt')
```

We can create a SPARC-X-API calculator using similar parameters. Note
that in real-space DFT, the parameter mesh spacing (`h`) controls the
convergence. To avoid large "egg-box effect" due to large mesh
spacing, we recommend to use a smaller `h` value. For demonstration
purpose a rather rough mesh spacing `h=0.25` (in Angstrom) and a 3x3x3
k-points are used.

```{code} python
import numpy as np
from ase.build import bulk
from ase.eos import calculate_eos

from sparc import SPARC


def main():
    # Al in conventional cell
    atoms = bulk("Al", cubic=True)
    calc = SPARC(h=0.25, kpts=(3, 3, 3), xc="pbe", directory="ex0-eos")
    vol = atoms.get_volume()
    atoms.calc = calc
    eos = calculate_eos(atoms, npoints=5, eps=0.05, trajectory="al-eos-sparc.traj")
    print("Original volume: Ang^3", vol)
    v, e, B = eos.fit()
    print("Fitted volume (Ang^3), energy (eV), modulus (eV/Ang^3)")
    print(v, e, B)
    a0 = v ** (1 / 3)
    print(f"Optimal cell length (cubic): {a0} Ang")
    atoms.set_cell([a0, a0, a0], scale_atoms=True)
    e0_sparc = atoms.get_potential_energy()
    print(f"Energy calculated by SPARC: {e0_sparc} eV")
    print(f"Energy diff {abs(e0_sparc - e)} eV")
    return


if __name__ == "__main__":
    main()
```

```{note}
This example uses file I/O mode for demonstration purpose only. Consider choosing the [socket mode](advanced_socket.md) if you need more flexibility.
```

## Geometric optimization in file-I/O and socket modes

When a DFT workflow requires multiple evaluations of single point
calculations, running SPARC-X-API in socket mode will usually have
advantage over the standard file-I/O mode. Here we use an example of
optimizing an ammonia molecule to show the difference.

First we construct a NH3 molecule in a box with Dirichlet boundary
conditions, and optimize it with SPARC's internal geometric optimization
(geopt) routine. We use a rough mesh spacing for demonstration purpose only.
```{code} python
import numpy as np
from ase.build import molecule
from ase.constraints import FixAtoms
from sparc import SPARC

nh3 = molecule("NH3", cell=(8, 8, 8), pbc=False)
# Fix the N center
nh3.constraints = [FixAtoms([0])]
nh3.rattle()


def optimize_sparc_internal():
    atoms = nh3.copy()
    calc = SPARC(
        h=0.25,
        kpts=(1, 1, 1),
        xc="pbe",
        convergence={"forces": 0.02},
        relax_flag=True,
        print_relaxout=True,
        relax_method="LBFGS",
        directory="ex1-sparc",
    )
    atoms.calc = calc
    e_fin = atoms.get_potential_energy()
    f_fin = atoms.get_forces()
    nsteps = len(calc.raw_results["geopt"])
    print("SPARC internal LBFGS:")
    print(f"Final energy: {e_fin} eV")
    print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    print(f"N steps: {nsteps}")
```

Alternatively we can use any of the `ase.optimize` optimizers to perform
the geometric relaxation, for example BFGS.
The following code uses the file I/O mode:
```{code} python
def optimize_ase_bfgs():
    atoms = nh3.copy()
    calc = SPARC(
        h=0.25,
        kpts=(1, 1, 1),
        xc="pbe",
        directory="ex1-ase"
    )
    atoms.calc = calc
    opt = BFGS(atoms)
    opt.run(fmax=0.02)
    e_fin = atoms.get_potential_energy()
    f_fin = atoms.get_forces()
    nsteps = opt.nsteps
    print("ASE LBFGS")
    print(f"Final energy: {e_fin} eV")
    print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    print(f"N steps: {nsteps}")
```

There are several drawbacks in file I/O mode with multiple single point
calculations:

1) The number of `.static` files can accumulate quickly (e.g. `.static_01`, `.static_02`, etc)
2) In each calculation the density and orbital are re-calculated
3) There are overhead when writing / loading files

You can overcome these using the socket mode, effectively by just
invoking the `use_socket=True` flag in above case.
```{note}
Make sure your SPARC binary is compiled with socket support. See [installation guide](installation.md) for more details.
```

```{code} python
def optimize_ase_bfgs_socket():
    atoms = nh3.copy()
    calc = SPARC(
        h=0.25, kpts=(1, 1, 1), xc="pbe", print_forces=True, directory="ex1-ase-socket",
        use_socket=True,
    )
    atoms.calc = calc
    with calc:
        opt = BFGS(atoms)
        opt.run(fmax=0.02)
    e_fin = atoms.get_potential_energy()
    f_fin = atoms.get_forces()
    nsteps = opt.nsteps
    print("ASE LBFGS")
    print(f"Final energy: {e_fin} eV")
    print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    print(f"N steps: {nsteps}")
```



## Training Machine Learned Force Fields using SPARC and the SPARC-X-API

The SPARC source code has the ability to train MLFF on-the-fly using
the SOAP descriptor and Bayesian Linear Regression. This feature can
be accessed via the API by passing the appropriate sparc flags to the
`SPARC` calculator object. An exhaustive list of MLFF parameters are
available in the [SPARC
documentation](https://github.com/SPARC-X/SPARC/tree/master/doc).

The primary flag of interest is the `MLFF_FLAG` which must be set to `1` to train a MLFF from scratch. This needs to be included in a parameter dictionary that will later be passed to a SPARC calculator instance:

```python
calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1,1,1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "ELEC_TEMP_TYPE": "fermi-dirac",
    "ELEC_TEMP": 116,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "MLFF_FLAG": 1,
    "MLFF_INITIAL_STEPS_TRAIN": 3,
}
```
This parameter dictionary primes an on-the-fly simulation. The `MLFF_INITIAL_STEPS_TRAIN` keyword specifies the number of reference structures that will be added to the training set before the model is first trained. References structures can be generated using the SPARC calculator in FileIO mode with C SPARC's internal relaxation or MD algorithms. Adding

```python
calc_params['RELAX_FLAG'] = 1
```

or

```python
calc_params['MD_FLAG'] = 1
```

to the parameter dictionary will train an ML model on-the-fly during relaxation or MD. Additional MD related parameters will be necessary such as the timestep and method. See the SPARC documentation for details. The simulations can be triggered by calling the familiar property functions on an `ASE` `Atoms` object. Here, we use a water molecule:

```python
from sparc.calculator import SPARC
from ase.build import molecule

water = molecule('H2O', vacuum=7)
water.pbc = [False,False,False]
water.calc = SPARC(**calc_params)
energy = water.get_potential_energy()
```

This triggers an on-the-fly relaxation or MD simulation. Additional details are available in the SPARC documentation regarding hyperparameters and recommended settings.

## Pairing SPARC and the SPARC-X-API with external algorithms via iPI socket communication

The iPI communication protocol embedded in SPARC and accessed via the SPARC-X-API allows users to treat SPARC as a DFT backend for many additional python libraries. For example, the reference structures for training MLFF on-the-fly can be generated via external algorithms such as `ASE` optimizers or MD engines. The DFT energies, forces, and stresses are computed by SPARC and parsed by the API to allow positions to be updated externally. The code from the previous section can be modified to use the `ASE` implementation of the BFGS algorithm:

```python
from sparc.calculator import SPARC
from ase.build import molecule
from ase.optimize import BFGS

water = molecule('H2O', vacuum=7)
water.pbc = [False,False,False]

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1,1,1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "ELEC_TEMP_TYPE": "fermi-dirac",
    "ELEC_TEMP": 116,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "MLFF_FLAG": 1,
    "MLFF_INITIAL_STEPS_TRAIN": 3,
}

with SPARC(use_socket=True, **calc_params) as calc:
    water.calc = calc
    dyn = BFGS(water, trajectory = 'water-bfgs-opt.traj')
    dyn.run(fmax=0.05)
```

A MLFF is still trained on-the-fly using SPARC DFT energies, forces, and stresses; but the structures are generated from the `ASE` optimization algorithm.

## Advanced Usage: Training MLFF on-the-fly during metadynamics simulation

The socket makes porting SPARC to external codes with python interfaces trivial. The PLUMED pacakge's existing `ASE` interface can be coupled with the SPARC-X-API to train MLFF on-the-fly during metadynamics simulations. The setup on the SPARC side remains unchanged from other socket examples:

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
