# sparc-dft-api
[![Package](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/package.svg)](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/package.svg)
[![Coverage](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/coverage.svg)](https://raw.githubusercontent.com/alchem0x2A/sparc-dft-api/badges/badges/coverage.svg)
[![Unit tests](https://github.com/alchem0x2A/sparc-dft-api/actions/workflows/installation_test.yml/badge.svg)](https://github.com/alchem0x2A/sparc-dft-api/actions/workflows/installation_test.yml)

`sparc-dft-api` is an [ASE](https://wiki.fysik.dtu.dk/ase/)-compatible python API for the density functional theory (DFT) code [SPARC](https://github.com/SPARC-X/SPARC). Starting v0.2, it provides the following functionalities:

1. ASE-compatible I/O format for SPARC files
2. JSON API associated with SPARC C-code for parameter validation and conversion
3. Fully functional calculator interface for SPARC

<!-- *TODO*:
- [ ] More advanced interface
- [ ] 
 -->
## Installation:

`sparc-dft-api` may be installed via any of the following approaches:

### 1. Via `anaconda` or `miniconda` (recommended)

Install `anaconda` or `miniconda` and create a working conda environment (see [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)). After than, use `conda install` to install the API:

```bash
conda install -c alchem0x2a sparc-dft-api
```
<!-- *TODO*:
- [] Change to conda-forge later
 -->

*NOTE: above conda code will change to conda-forge channel in official release*


You may also want to install our pre-compiled SPARC binaries and the SPMS pseudopotentials within the environment (x86-64 Linux only)

```bash
conda install -c alchem0x2a sparc
```

*NOTE: above conda code will change to conda-forge channel in official release*


<!-- *TODO*:
- [ ] Push the conda-package to channel
- [ ] Make sure the SPARC_PP_PATH env variables are automatically set
- [ ] Mantain a conda-forge release?
 -->
 
### 2. Install stable version from [PyPI]()
```bash
python -m pip instal sparc-dft-api>=2.0
```

*NOTE: need to update compiled wheel in pypi*

<!-- WIP
- [ ] Push to pypi? 
-->

### 3. Install latest develop version from GitHub

```bash
python -m pip install git+https://github.com/SPARC-X/sparc-dft-api
```

You can download the latest SPMS pseudopotentials after installation (optional):
```bash
python -m sparc.download_data
```

Please following the SPARC [manual](https://github.com/SPARC-X/SPARC) for compilation and installation of the SPARC DFT code itself in this case.

<!-- *TODO*
- [ ] Make Pypi
- [ ] Make wheel available
 -->

## Setting up the environment
By design, `sparc-dft-api` >= v0.2 will automate the discovery for 
pseudopotential files, JSON API and SPARC binary. 
But you can also have fine control over how they can be setup:

### A) Pseudopotential files
Pseudopotential files (in `Abinit` psp8 format) are looked for in the following
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

*Note*: the `-name PREFIX` part can be omitted in `sparc-dft-api` >= v0.2, since it will autocomplete the command using the `SPARC.label` property.

## Post-installation check

We recommend the users to run a simple check after installation and setup:
```bash
python -m sparc.quicktest
```

A successful setup would show the following blocks at the end of the output

<img width="500" alt="image" src="https://github.com/alchem0x2A/sparc-dft-api/assets/6829706/95cb712e-4c77-4b14-8130-4961e3c50278">



If you encounter any issues, please refer to the [Trouble Shooting]() section.

<!-- *TODO*
- [ ] Test ase io format support
- [ ] Test calculator command
- [ ] Test API version
- [ ] Complete trouble shooting
 -->


## `sparc-dft-api`: Basic usages
### 1. Read / write SPARC files

Unlike other DFT codes where the ASE format corresponds to a single file,
`sparc-dft-api` provides I/O support for the whole calculation directory, a.k.a "SPARC bundle".``
`sparc-dft-api` allows automatic discovery of this file format (`sparc-bundle`) in ASE:

- Read from a sparc bundle

```python
import sparc
from ase.io import read, write

atoms = read("test.sparc", format="sparc", image=-1)
```

- Write a minimal sparc bundle from atoms

```python
import sparc
from ase.io import read, write
from ase.build import Bulk
atoms = Bulk("Al") * [4, 4, 4]
atoms.write("test.sparc", format="sparc")
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

This example sets up a calculation for H2 atoms in a 10 x 10 x 10 Å$^3$ PBC cell with PBE 
exchange correlation function, and a grid spacing (`h`) of 0.25 Å. Note by calling `atoms.get_forces`,
the calculator will automatically sets the flags for printing the forces.

2. Geometric optimization (using SPARC's internal routines)
```python
from sparc.calculator import SPARC
from ase.build import bulk
atoms = bulk("Al", cubic=True)
atoms.rattle()
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


## Major changes from `sparc-dft-api` [v0.1](https://github.com/SPARC-X/sparc-dft-api/tree/eac557f214b402122a506f88f38c7a8767283503)

`sparc-dft-api` has been heavily refactored in v0.2. If you're using legacy Python codes
that are written under v0.1 API, there are a few major changes that require your attention:

1. Support for single `.ion` file format is deprecated. Instead, `v0.2` API treats the whole SPARC directory as a bundle format. Please use `read_sparc` and `write_sparc` methods for basic file I/O instead.
Nevertheless, reading calculation results generated by a v0.1 API code will not be affected.

2. v0.2 API uses a different mapping scheme for the sorting of ASE atoms objects (similar to `Vasp`), add a comment section in `.ion` file similar to follows:
```python
# ASE-SORT:
# 3 2 1 0
# END ASE-SORT
```
which maps atoms 3, 2, 1, 0 from the SPARC .ion file order to atoms 0, 1, 2, 3 in ASE order. This is useful for systems that are constructed by ASE's `add_adsorbate` method.

3. v0.2 API accepts all SPARC internal parameters (i.e. **CAPITALIZED**) in *atomic units* for consistency reason. 
However, we also keep a list of "special input params" that are conventionally used in other ASE calculators, that use Å / eV / GPa / fs unit system.

4. Defining `LATVEC`, `LATVEC_SCALE`, or `CELL` via the calculator parameters is no longer encouraged. Instead, all structure changes should be made to the `Atoms` object.

For more discussion please see [Advanced Topic] section.

Below are a list of v0.1 method of the `SPARC` calculator and their current status in v0.2 API. 
`calc` is an instance of `sparc.SPARC`.

|  old methods           | status in v0.2 API |     alternatives                   |
|------------------------|--------------------|------------------------------------|
| `interpret_grid_input` | deprecated         | `calc.set(fd_grid=[20, 20, 20])`   |
| `interpret_kpoint_input` | deprecated         | `calc.set(kpts=[3, 3, 3])`   |
| `interpret_downsampling_input` | deprecated   | Manual setting not recommended |
| `interpret_kpoint_shift` | deprecated | `calc.set(kpoint_shift=[0, 0, 0])` |
| `get_pseudopotential_directory` | deprecated | `calc.psp_dir` |
| `get_nstates`          | maintained     |            |
| `setup_parallel_env` | deprecated | Manual set |
| `generate_command` | deprecated | `calc._make_command()` |
| `estimate_memory` | maintained | |
| `get_scf_steps` | maintained | |
| `get_geometric_steps` | deprecated | `calc.get_number_of_ionic_steps()`|
| `get_runtime` | maintained | |
| `get_fermi_level` | maintained | |
| `concatinate_output` | deprecated | Use `sparc.SparcBundle` instead |
| `read_line` | deprecated | Use `sparc.SparcBundle` instead |
| `parse_output` | deprecated | `calc.read_results()` |
| `parse_relax` | deprecated | `calc.read_results()` |
| `parse_md` | deprecated | `calc.read_results()` |
| `parse_input_args` | deprecated | `calc.set(**kwargs)` |
| `recover_index_order_from_ion_file` | deprecated | Use `calc.sort` and `calc.resort` |
| `atoms_dict` | deprecated | Use third party library like `bson` |
| `dict_atoms` | deprecated | Use third party library like `bson` |

## Advanced Topics
### Behind the bundle file format

Instead of parsing individual `.ion` and `.inpt` files, 
the bundle format (recognized by ASE by `format="sparc"`) will 
gather information from all files and check if atomic information
and calculation results can be retrieved. 
The central piece for handling the bundle format is 
`sparc.io.SpardBundle` class. You can use it to parse an existing bundle

```python
from sparc.io import SparcBundle
bundle = SparcBundle("path/to/your-calc.sparc", mode="r")
images = bundle.convert_to_ase(index=":")
```

or write an `Atoms` object to a bundle with a minimal set of `.ion` and `.inpt`
files.

```python
from sparc.io import SparcBundle
from ase.io import read
atoms = read("some-externalfile.xyz")
bundle = SparcBundle("path/to/your-calc.sparc", mode="w")
bundle._write_ion_and_inpt(atoms, label="SPARC")
```

For each individual SPARC file (e.g. `.ion`, `.inpt`, `.static`, `.geopt`, `.aimd`), 
file-specific parsers are in `sparc.sparc_parsers.<format>` files. 
Each `_read_<format>` method will return the structured raw-data dictionary of the files.
Similarly, `_write_<format>` takes the structured dictionary as input and write the file
using only relevant data.

### Behind the JSON API

The JSON API are directly parsed from the `SPARC` documentation [LaTeX files](https://github.com/SPARC-X/SPARC/tree/master/doc/.LaTeX). 
The JSON API file (`sparc/sparc_json_api/parameters.json`) distributed by `sparc-dft-api` is generated by:

```bash
git clone https://github.com/SPARC-X/SPARC.git
python -m sparc.docparser SPARC/doc/.LaTeX
```

You can use `sparc.api.SparcAPI` together with the generated JSON API schema:
1. Load the latest SPARC api
```python
from sparc.api import SparcAPI
sis = SparcAPI()
print(sis.sparc_version)
```


2. Check if a variable is available
```python
from sparc.api import SparcAPI
sis = SparcAPI()
# A typo will be detected (actual parameter is CALC_PRESS)
assert "CALC_PRESSURE" not in sis.parameters
```

3. Convert string <--> valid value
```python
from sparc.api import SparcAPI
sis = SparcAPI()
latvec = sis.convert_string_to_value("LATVEC", "1.0 0 0\n 0 1.0 0\n 0 0 1.0")
latvec_string = sis.convert_value_to_string("LATVEC", latvec)
```

4. Provide help info for a given parameter
```python
from sparc.api import SparcAPI
sis = SparcAPI()
```

### Retriving parameters from old SPARC calculations

`sparc.SPARC` calculator supports the `restart` mode which will reconstruct all 
parameters, psp files and atomic information from an existing SPARC calculation and 
rerun them. 

```python
from sparc import SPARC
calc = SPARC(restart=True, directory="old-calc.sparc")
old_atoms = calc.atoms.copy()
# Redo the calculation with updated parameters
old_atoms.rattle()
calc.set(h=0.2, directory="new-calc.sparc")
old_atoms.calc = calc
old_atoms.get_potential_energy()
```


### Special inputs for `sparc.SPARC` calculator

The following input parameters have special meaning in `sparc.SPARC` calculator,
they are consistent with definitions in other ASE calculators and uses Å / eV / GPa / fs 
unit system:

| parameter name | meaning                         | example        | equivalent `SPARC` input |
|----------------|---------------------------------|----------------|--------------------------|
| `xc`           | Exchange-correlation functional | `xc=pbe` | `EXCHANGE_CORRELATION: GGA_PBE` |
| `h`            | Real grid spacing    (Å)           | `h=0.2`        | `FD_GRID: Nx Ny Nz`  (calculated values)         |
| `gpts`         | Explicit grid points |   `gpts=[10, 10, 10]` | `FD_GRID: 10 10 10` |
| `kpts`         | Kpoint mesh          |   `kpts=[3, 3, 3]`    | `KPOINT_GRID: 3 3 3` |
| `convergence`  | Dict of convergence criteria (see below) |  | |
|                | `energy`  eV/atom         | `convergence={"energy": 1e-4}` | `SCF_ENERGY_ACC: 3e-6` |
|                | `forces`  eV/Å            | `convergence={"forces": 1e-2}` | `TOL_RELAX: 2e-4` |
|                | `density` e/atom          | `convergence={`density`: 1e-6}`| `TOL_PSEUDOCHARGE: 1e-6` |

*WIP*: support more advanced settings like D3, HSE and DFT+U in `xc` settings
