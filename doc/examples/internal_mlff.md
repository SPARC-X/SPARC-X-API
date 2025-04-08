# Using Machine learning force fields (MLFF) in SPARC-X-API

## Training Machine Learned Force Fields using SPARC and the SPARC-X-API

The SPARC source code has the ability to train MLFF on-the-fly using
the SOAP descriptor and Bayesian Linear Regression. This feature can
be accessed via the API by passing the appropriate sparc flags to the
`SPARC` calculator object. An exhaustive list of MLFF parameters are
available in the [SPARC
documentation](https://github.com/SPARC-X/SPARC/tree/master/doc).

The full python script for training MLFF using SPARC calculated geometric relaxation data can be downloaded [here](https://github.com/SPARC-X/SPARC-X-API/blob/master/examples/FileIO/relax/relax_coords/mlff/run.py).

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

to the parameter dictionary will train an ML model on-the-fly during
relaxation or MD. Additional MD related parameters will be necessary
such as the timestep and method. See the SPARC documentation for
details. The simulations can be triggered by calling the familiar
property functions on an `ASE` `Atoms` object. Here, we use a water
molecule:

```python
from sparc.calculator import SPARC
from ase.build import molecule

water = molecule('H2O', vacuum=7)
water.pbc = [False,False,False]
water.calc = SPARC(**calc_params)
energy = water.get_potential_energy()
```

This triggers an on-the-fly relaxation or MD simulation. Additional
details are available in the SPARC documentation regarding
hyperparameters and recommended settings.

## Pairing SPARC and the SPARC-X-API with external algorithms via iPI socket communication

The iPI communication protocol embedded in SPARC and accessed via the
SPARC-X-API allows users to treat SPARC as a DFT backend for many
additional python libraries. For example, the reference structures for
training MLFF on-the-fly can be generated via external algorithms such
as `ASE` optimizers or MD engines. The DFT energies, forces, and
stresses are computed by SPARC and parsed by the API to allow
positions to be updated externally. The code from the previous section
can be modified to use the `ASE` implementation of the BFGS algorithm:

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

The full script for socket-enabled MLFF training and optimization can be downloaded [here](https://github.com/SPARC-X/SPARC-X-API/blob/master/examples/socket/relax/mlff/coords/run.py).
