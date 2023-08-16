# `sparc-dft-api`: A Python API for the SPARC DFT Code
[![Package](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/package.svg)](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/package.svg)
[![Coverage](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/coverage.svg)](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/coverage.svg)
[![Unit tests](https://github.com/alchem0x2A/sparc-dft-api/actions/workflows/installation_test.yml/badge.svg)](https://github.com/alchem0x2A/sparc-dft-api/actions/workflows/installation_test.yml)

`sparc-dft-api` is an [ASE](https://wiki.fysik.dtu.dk/ase/)-compatible Python API for the density functional theory (DFT) code [SPARC](https://github.com/SPARC-X/SPARC). It offers:

1. ASE-compatible I/O format for SPARC files
2. A JSON API interfacing with SPARC's C-code for parameter validation and conversion
3. A comprehensive calculator interface for SPARC.


## Installation:

The Python API may be installed via either of the following approaches:

### 1. Via `anaconda` or `miniconda` (recommended)

Set up a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) and install the Python API,
which includes the pseudopotential files:

```bash
# Change 'sparc-env' to your desired name if needed
conda create -n sparc-env
conda activate sparc-env
conda install -c alchem0x2a sparc-dft-api
```

On Linux platforms (x86_64, aarch64), you can also install the
precompiled `sparc` DFT binaries alongside the API:

```bash
conda install -c alchem0x2a sparc
conda activate sparc-env   # Re-activate to have the env variables effective
```

*Note: Packaging of sparc-dft-api on conda-forge is in progress.*

### 2. Manual installation from source with `pip`


```bash
python -m pip install git+https://github.com/SPARC-X/sparc-dft-api
```

Optionally, you can download the latest SPMS pseudopotentials
post-installation, if you don't have them already:

```bash python -m sparc.download_data ```

This command unpacks the pseudopotential files into
`<python-lib-root>/site-packages/sparc/psp`.

To utilize the API for initiating SPARC calculations, please
following the [SPARC manual](https://github.com/SPARC-X/SPARC) for
compilation and installation of the SPARC DFT code itself.

## Post-installation check

We recommend the users to run a simple test after installation and setup:

```bash
python -m sparc.quicktest
```

A proper setup will display the following sections at the output's conclusion:

<img width="500" alt="image" src="https://github.com/alchem0x2A/sparc-dft-api/assets/6829706/95cb712e-4c77-4b14-8130-4961e3c50278">

For using the API to interpret SPARC input and output files, it's
essential that the "Import" and "JSON API" tests are successful. For
calculations to be viable, all tests must pass.

Please refer to the [Setting Up the
Environment](#setting-up-the-environment) or guidance on correctly
configuring the environment variables. This ensures the Python API can
locate the SPARC setups.  If you run into any problems, consult our
[Trouble Shooting](doc/troubleshooting.md).

## Setting up the environment
By design, `sparc-dft-api` will automate the discovery for
pseudopotential files, JSON API and SPARC binary.
But you can also have fine control over how they can be setup:

### A) Pseudopotential files
Pseudopotential files (in `Abinit` psp8 format) are searched in the following
order:
1) `psp_dir` argument passed to the `sparc.SPARC` calculator
2) Environmental variables `$SPARC_PSP_PATH` or `$SPARC_PP_PATH`
3) `psp8` files bundled with the sparc-dft-api installation (must be downloaded via `python -m sparc.download_data`)

You can  set the `$SPARC_PSP_PATH` variable like follows:
```bash
export SPARC_PSP_PATH="/path/to/your/psp8/directory"
```

To get the location of default psp8 files in option 3), run the following code:
```bash
python -c "from sparc.common import psp_dir; print(psp_dir)"
```

### B) JSON API file
Currently the calculator's API validator is bundled with sparc-dft-api.
In future releases it will be possible to dynamically load a JSON API by
matching the SPARC binary version. You can take a look at the api file at
`sparc.sparc_json_api.default_json_api`, which should contain entries like:

```json
"FD_GRID": {
   "symbol": "FD_GRID",
   "label": "FD_GRID",
   "type": "integer array",
   "default": null,
   "unit": "No unit",
   "example": "FD_GRID: 26 26 30",
   "description": "A set of three whitespace delimited values specifying the number of finite-difference intervals in the lattice vector (LATVEC) directions, respectively.",
   "remark": "The convergence of results with respect to spatial discretization needs to be verified. ECUT, MESH_SPACING, FD_GRID cannot be specified simultaneously.",
   "allow_bool_input": false,
   "default_remark": "None",
   "description_raw": "A set of three whitespace delimited values specifying the number of finite-difference intervals in the lattice vector (\\hyperlink{LATVEC}{\\texttt{LATVEC}}) directions, respectively.",
   "remark_raw": "The convergence of results with respect to spatial discretization needs to be verified. \\hyperlink{ECUT}{\\texttt{ECUT}}, \\hyperlink{MESH_SPACING}{\\texttt{MESH\\_SPACING}}, \\hyperlink{FD_GRID}{\\texttt{FD\\_GRID}} cannot be specified simultaneously.",
   "category": "system"
  },
```

### C) `SPARC` command

The command used for running SPARC calculations are detected in the following order:
1) Passing `command` argument to `sparc.SPARC` calculator
2) Setting `$ASE_SPARC_COMMAND` variable
3) If None of the above exists, look for a SPARC binary under current `$PATH` and combine with the suitable `mpi` command prefix (auto-detected).

Example:
```bash
export ASE_SPARC_COMMAND="mpirun -n 8 /path/to/sparc -name PREFIX"
```

*Note*: the `-name PREFIX` part can be omitted in `sparc-dft-api` > v1.0.0, since it will autocomplete the command using the `SPARC.label` property.



## `sparc-dft-api`: Basic usages
### 1. Read / write SPARC files

Unlike other DFT codes where the ASE format corresponds to a single file,
`sparc-dft-api` provides I/O support for the whole calculation directory, a.k.a "SPARC bundle".
`sparc-dft-api` allows automatic discovery of this file format (`"sparc"`) in ASE:

- Read from a sparc bundle

```python
import sparc
from ase.io import read, write

atoms = read("test.sparc", index=-1)
```

- Write a minimal sparc bundle from atoms

```python
import sparc
from ase.io import read, write
from ase.build import Bulk
atoms = Bulk("Al") * [4, 4, 4]
atoms.write("test.sparc")
```

For more details about the bundle IO format, please see [Advanced Topics]() section

### 2. JSON API for SPARC calculator

ASE-based DFT calculator interfaces often face the challenge of version inconsistencies
between low-level codes (Fortran/C/C++) and outdated upper-level APIs (Python) regarding parameter sets and default values.
In `sparc-dft-api`, DFT parameters are meticulously managed through a JSON API, translated from SPARC's LaTeX documentation.
Each release of `sparc-dft-api` is linked with a specific version of the SPARC source code, ensuring compatibility and consistency with the default parameter set.

Automatic conversion and validation between plain-text SPARC input files and sparc-dft-api are
handled using the `sparc.api.SparcAPI` class. If you want more fine control, please refer to [Advanced Topics]() section.


### 3. Calculator interface

`sparc-dft-api` provides a calculator interface that should be familiar for users with experience on any other
ASE calculators (`Vasp`, `QuantumEspresso`, `GPAW` etc), as shown in the following examples:

1. Single point calculation
```python
from sparc.calculator import SPARC
from ase.build import molecule
atoms = molecule("H2", cell=(10, 10, 10), pbc=True, directory="run_sp")
atoms.calc = SPARC(h=0.25)
atoms.get_potential_energy()
atoms.get_forces()
```

This example sets up a calculation for H2 atoms in a 10 Å x 10 Å x 10 Å PBC cell with PBE
exchange correlation function, and a grid spacing (`h`) of 0.25 Å. Note by calling `atoms.get_forces`,
the calculator will automatically sets the flags for printing the forces.

2. Geometric optimization (using SPARC's internal routines)
```python
from sparc.calculator import SPARC
from ase.build import bulk
atoms = bulk("Al", cubic=True)
atoms.rattle(0.05)
atoms.calc = SPARC(h=0.25, kpts=(3, 3, 3), relax_flag=True, directory="run_opt")
atoms.get_potential_energy()
atoms.get_forces()
```

This example sets up a calculation for a rattled Aluminum primitive unit cell, calculate with PBE
functional, grid spacing of 0.25 Å, and 3 x 3 x 3 K-point grid.
Optimization of ionic positions
is handled with SPARC's default LBFGS routine.

3. AIMD in SPARC

```python
from sparc.calculator import SPARC
from ase.build import bulk
md_params = dict(md_flag=True, ion_temp=800, md_method="NVE", md_timestep=0.6, md_nstep=5)
atoms = bulk("Al") * (3, 3, 3)
atoms.rattle()
atoms.calc = SPARC(h=0.25, kpts=(1, 1, 1), directory="run_aimd", **md_params)
atoms.get_potential_energy()
```

This example runs a short NVE MD simulation (5 steps) at 800 K for 27 Al atoms.
If you want to extract more information about the MD simulation steps, take a look at `SPARC.raw_results`.

4. Geometric optimization using ASE's optimizers

The power of `sparc-dft-api` is to combine single point `SPARC` calculations with advanced ASE optimizers, such as BFGS, FIRE or GPMin. Example 2 can be re-written as:

```python
from sparc.calculator import SPARC
from ase.build import bulk
from ase.optimize import LBFGS
atoms = bulk("Al", cubic=True)
atoms.rattle(0.05)
atoms.calc = SPARC(h=0.25, kpts=(3, 3, 3), directory="run_opt_ase")
opt = LBFGS(atoms, alpha=90)
opt.run(fmax=0.02)
```

### 4. Commandline tools

`sparc-dft-api` provides a simple command wrapper `sparc-ase` to add
support of SPARC file formats to the `ase` cli tools. Simple replace
`ase [subcommand] [args]` with `sparc-ase [subcommand] [args]` to
access your SPARC bundle files as you would use for other file formats.
As an example, use `sparc-ase gui path/to/your/bundle.sparc` for the visualization of
atomistic structures. Depending on the file contents, either single atoms or multiple
images will be displayed.

Below is a screenshot showing the usage of `sparc-ase gui` to visualize a
short [MD trajectory](tests/outputs/NH3_sort_lbfgs_opt.sparc).

<img width="1200" alt="image" src="https://github.com/alchem0x2A/sparc-dft-api/assets/6829706/e72329ff-7194-4819-94f8-486ef2218844">

## Troubleshooting
Please refer to the [troubleshooting](doc/troubleshooting.md) guidelines

## Advanced topics
A detailed description about how the API works can be found [here](doc/advanced_topics.md)

## API changes
The API changes compared to v0.1 are summarized [here](doc/changes_v0.1.md)

## How to contribute
Please refer to our [guidelines for contributors](doc/contribution_guideline.md)
