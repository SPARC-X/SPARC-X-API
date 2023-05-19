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
### 0. Prerequisites

By design, sparc-dft-api >= v0.2 will automate the discovery for 
pseudopotential files, JSON API and SPARC binary. 
But you can also have fine control over how they can be setup:

### A) Pseudopotential files
Pseudopotential files (in *Abinit*'s psp8 format) are looked for in the following
order:
1) `psp_dir` argument passed to the `sparc.SPARC` calculator
2) Environmental variables `$SPARC_PSP_PATH` or `$SPARC_PP_PATH` 
3) `psp` directory bundled with the sparc-dft-api installation (must be downloaded via `python -m sparc.download_data`)

### B) JSON API file
Currently the calculator's API validator is bundled with sparc-dft-api. 
In future releases it will be possible to dynamically load a JSON API by 
matching the SPARC binary version.

### C) `SPARC` command

The command used for running SPARC calculations are detected in the following order:
1) Passing `command` argument to `sparc.SPARC` calculator
2) Setting `$ASE_SPARC_COMMAND` variable
3) If None of the above exists, look for a SPARC binary under current `$PATH` and combine with the suitable `mpi` command prefix (auto-detected).

Example:
```bash
export ASE_SPARC_COMMAND="mpirun -n 8 /path/to/sparc -name NAME"
```

*Note*: the `-name NAME` part can be omitted. In that case, the calculator will 
autocomplete the command using its `label` property.


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
