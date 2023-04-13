# sparc-dft-api

sparc-dft-api is an ASE based python wrapper for the density functional theory (DFT) code SPARC. This wrapper requires <ins>Python3</ins>, and is *currently in active development and alpha testing*, so please use it with caution and report any errors in the "Issues" and/or submit pull requests with patches to fix issues as needed.

# Work-in-progress

The current sparc python API is under heavy re-construction. The below are the on-going API changes


# Allowable Arguments

`sparc` Python-API are directly parsed from the `SPARC` c-code's documentation LaTeX files. 
We provide a document parser to extract the allowed parameters in the current version of `SPARC`.

```
python -m sparc.docparser <root-to-sparc-cpde>/doc/.LaTeX
```

which generates a json file `parameters.json` for the current API. 
We also keep track of the latest API under `sparc/sparc_json_api/parameters.json`.

We provide a API object `sparc.inputs.SparcInputs` for loading and validating parameters.

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
~~~

