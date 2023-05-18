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

### Post-installation check

We recommend the users to run a simple check after installation:
```bash
python -m sparc.quicktest
```

A successful setup would have the following output
**TODO: add the image placeholder**

If you encounter any issues, please refer to the [Trouble Shooting] section.

<!-- *TODO*
- [ ] Test ase io format support
- [ ] Test calculator command
- [ ] Test API version
- [ ] Complete trouble shooting
 -->


## `sparc-dft-api`: Basic usages
### 0. Environmental variables
By design, sparc-dft-api >= v0.2 will automate the discovery for 
pseudopotential files, JSON API and SPARC binary. You can control them 
by setting the environmental varialbles:
### A) Pseudopotential files
Pseudopotential files (in *Abinit*'s psp8 format) are looked for in the following
order:
1) `psp_dir` argument passed to the `sparc.SPARC` calculator
2) Environmental variables `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` 
3) `psp` directory bundled with the sparc-dft-api installation (must be downloaded via `python -m sparc.download_data`)

### B) JSON API file
**TODO**

### C) `SPARC` command
**TODO**


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

- Write a minimal sparc bundle from atoms / images

```python
import sparc
from ase.io import read, write
from ase.build import Bulk
atoms = Bulk("Al") * [4, 4, 4]
atoms.write("test.sparc", format="sparc")
```

For more details about the bundle IO format, please see [Behind the scene] section

### 2. JSON API for SPARC calculator

A common problem in quantum mechanical codes is that the low-level codes (written in Fortran/C/C++) may 
introduce changes from version to version, while the upper-level API (ASE calculator interface) 
are often out-dated regarding the parameter sets, default values etc.

DFT parameters in `sparc-dft-api` are *first-world members*, by using a JSON API translated from the LaTeX
documentation of SPARC. Each release of `sparc-dft-api` will have the default parameter set linked with the same 
version of SPARC executable, while providing compatibility with other versions.

We provide a API object `sparc.inputs.SparcInputs` for loading and validating parameters, as well as generating 
user-friendly help information. Most most 
users, the API validation / conversion is done automatically. 
If you want more control over the JSON API, please refer to the [JSON API] section.

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
atoms = bulk("Al")
atoms.rattle()
atoms.calc = SPARC(h=0.25, kpts=(3, 3, 3), relax_flag=True, directory="run_opt")
atoms.get_potential_energy()
atoms.get_forces()
```

This example sets up a calculation for a rattled Aluminum primitive unit cell, calculate with PBE 
functional, grid spacing of 0.25 Å, and 3 x 3 x 3 K-point grid. Optimization of ionic positions 
is handled with SPARC's internal LBFGS routine.

3. AIMD in SPARC
*WIP*

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




## Major changes from `sparc-dft-api` [v1.0]

The python API has been heavily re-formatted v2.0. If you're using legacy python codes
that are written under v1.0 API, there are a few major changes that require your attention:

1. v2.0 API no longer provides a single I/O format for SPARC `.ion` file due to the reason explained in 
   [I/O] section. The lattice information in `.ion` files generated by legacy API will remain
   harmless but not parsed.
2. v2.0 API uses a different mapping scheme for the sorting of ASE atoms objects (similar to `Vasp`).
3. v2.0 API keeps all SPARC internal parameters (i.e. those can be **CAPITALIZED**) in atomic units, will all other
   ASE-specific units are in Å / eV / fs. *TODO: discuss with group for opinion*
4. v2.0 API is more flexible treating the `ASE_SPARC_COMMAND` environmental variable. While the same command v1.0 uses
   should still work, there is no need to specify `-name PREFIX` in the command.



<!-- ###### *Behind the scene*

`import sparc` will creates hooks in ASE's `ioformat` for `sparc.sparc_io_bundle.read_sparc` 
and `sparc.sparc_io_bundle.write_sparc` methods, allowing automatic format discovery.



## Behind the bundle file format

Instead of parsing individual .ion and .inpt files, 
the bundle format will gather information from all files and check if atomic information can be
retrieved.

New file-specific parsers will exist in `sparc.sparc_parsers.<format>` files. 
Each `_read_<format>` method will return the structured dictionary of the files.

Similarly, `_write_<format>` takes the structured dictionary as input and write the file
using only relevant data.

Implementations:
- [x] ion
- [x] inpt
- [ ] output
- [ ] geopt
- [ ] aimd
- [ ] multiple occurance
- [ ] Other files? Eigen? Grid results? 
  
## Allowable Arguments

`sparc` Python-API are directly parsed from the `SPARC` c-code's documentation LaTeX files. 
We provide a document parser to extract the allowed parameters in the current version of `SPARC`.

```
python -m sparc.docparser <root-to-sparc-cpde>/doc/.LaTeX
```

which generates a json file `parameters.json` for the current API. 
We also keep track of the latest API under `sparc/sparc_json_api/parameters.json`.



Use it like this:
1. Load the latest SPARC api
```python
from sparc.inputs import SparcInputs
sis = SparcInputs()
print(sis.sparc_version)
```

Output `2023.04.11`

2. Check if a variable is available
```python
from sparc.inputs import SparcInputs
sis = SparcInputs()
# A typo in 
assert "CALC_PRESSURE" not in sis.parameters
```

3. Convert string --> value and vice versa
```python
from sparc.inputs import SparcInputs
sis = SparcInputs()
# A typo in 
latvec = sis.convert_string_to_value("LATVEC", "1.0 0 0\n 0 1.0 0\n 0 0 1.0")
latvec_string = sis.convert_value_to_string("LATVEC", latvec)
print(latvec)
print(latvec_string)
```

4. Provide help info for a given parameter
```python
from sparc.inputs import SparcInputs
sis = SparcInputs()
print(sis.help_info("LATVEC"))
```

Output:
```
symbol: LATVEC
category: system
type: double array
unit: No unit
default: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
example: LATVEC: 0.5 0.5 0.0
0.0 0.5 0.5
0.5 0.0 0.5
description: A set of three vectors in row major order specifying the lattice vectors of the simulation domain (CELL).
remark: 
allow_bool_input: False
```

- [ ] Track the SPARC api for different versions 
- [ ] Allow `SparcInputs` loading different versions
- [ ] Allow calculator match the API version when loading the sparc binary
 



# OLD API

The old api specifications are from Ben Comer @GaTech.


## Installation:
The following command should install the package:

~~~
pip install git+https://github.com/SPARC-X/sparc-dft-api
~~~

prerequisites:

python packages:

`scipy`

`numpy`

`ase`

`spglib`

software:

`SPARC`

(https://github.com/SPARC-X/SPARC)

Currently installation is barebones, simply clone this repository (`git clone https://github.com/SPARC-X/sparc-dft-api.git`) and add the folder to your pythonpath (`export PYTHONPATH=[were you cloned the repository]:$PYTHONPATH`)


# Usage

To enable the SPARC calculator to work you need to tell it what command to use to run sparc and where your pseudopotentials are. These should be stored in the environment variables `$ASE_SPARC_COMMAND` and `$SPARC_PSP_PATH` respectively. e.g.:

`export ASE_SPARC_COMMAND=[command to run SPARC]`

`export SPARC_PSP_PATH=[location of pseudopotentials]`

The `ASE_SPARC_COMMAND` variable is somewhat tricky to set up. All the information about parallelization should be included where appropriate in the command (i.e. `mpirun -np X`.) SPARC takes in a flag `-name` followed by the prefix of the input file names. This wrapper deals with the prefixes completely internally, but the word `PREFIX` should be put in after the `-name` flag for SPARC. Here are two examples command if the `sparc` executable is in your `$PATH` environment variable:

`export ASE_SPARC_COMMAND="mpirun -np $PBS_NP sparc -name PREFIX"`


`export ASE_SPARC_COMMAND="mpirun -rmk pbs sparc -name PREFIX"`

note that the exact word "PREFIX" should be contained in this command.

`SPARC_PSP_PATH` is a path that points to the location of your pseudopotential files. These must have the naming convention: `[element name].pot` to be detected. Currently SPARC only uses the psp8 format. Alternatively, rather than setting this environment variable you may pass the information in with the `pseudo_dir` argument.

Once those environment varibles are in place, the SPARC ASE calculator works like any other ASE calculator. It must be imported, instantiated, and called. There are examples of common operations at the bottom of this README.

## Allowable Arguments

The SPARC calculator behaves similarly to other ASE calculators, taking an atoms object to define the system being calculated and a `label` argument for the filename prefixes. The arguments used to control the SPARC software are identical to those of the SPARC software program. Docuymentation can be found at the link below:

https://github.com/SPARC-X/SPARC/blob/master/doc/Manual.pdf

### Mesh Spacing
A mesh spacing or finite difference grid must be defined for the calculator to work. This can be done in one of three ways:
1. by inputting the MESH_SPACING argument. This will use SPARC's internal mesh spacing to generate a grid
2. by using the `h` argument, this  will default to converting to MESH_SPACING and utilizing this arguement in SPARC
3. Inputting the `fd_grid` arugment which explicitly defines the finite difference grid. For example `fd_grid=[25,25,25]`

## Examples


### Single Point Calculation
~~~
#get sparc calculator
from sparc.sparc_core import SPARC
calc = SPARC(h=0.2) # a grid spacing grid must be entered.

#make atoms
from ase.build import bulk
atoms = bulk('Si',cubic=True)
atoms.set_pbc([True] * 3)
atoms.set_calculator(calc)
atoms.get_potential_energy()
~~~


### Relaxation
~~~
from sparc.sparc_core import SPARC
calc = SPARC(h=0.2, RELAX_FLAG=1) # a grid spacing grid must be entered.

#make atoms
from ase.build import molecule
atoms = molecule('H2')
atoms.set_cell([6,6,6])
atoms.center()
atoms.set_pbc([False] * 3 )
atoms.set_calculator(calc)
atoms.get_potential_energy()
~~~

### Writing Input Files
~~~
from sparc.sparc_core import SPARC
calc = SPARC(h=0.2, RELAX_FLAG=1) # a grid spacing grid must be entered.

#make atoms
from ase.build import molecule
atoms = molecule('H2')
atoms.set_cell([6,6,6])
atoms.center()
atoms.set_pbc([False] * 3 )
atoms.set_calculator(calc)
calc.write_input()
~~~ -->

