# sparc-dft-api

sparc-dft-api is an ASE based python wrapper for the density functional theory (DFT) code SPARC. This wrapper is currently in alpha, so it's performance is not guaranteed.

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


## Usage

To enable the SPARC calculator to work you need to tell it what command to use to run sparc and where your pseudopotentials are. These should be stored in the environment variables `$ASE_SPARC_COMMAND` and `$SPARC_PSP_PATH` respectively. e.g.:

`export ASE_SPARC_COMMAND=[command to run SPARC]`

`export SPARC_PSP_PATH=[location of pseudopotentials]`

The `ASE_SPARC_COMMAND` variable is somewhat tricky to set up. All the information about parallelization should be included where appropriate in the command (i.e. `mpirun -np X`.) SPARC takes in a flag `-name` followed by the prefix of the input file names. This wrapper deals with the prefixes completely internally, but the word `PREFIX` should be put in after the `-name` flag for SPARC. Here is an example command if the `sparc` executable is in your `$PATH` environment variable:

`mpirun -np $PBS_NP sparc -name PREFIX`

Once those environment varibles are in place, the SPARC ASE calculator works like any other ASE calculator. It must be imported, instantiated, and called here is some example code for calculating bulk Si:

`SPARC_PSP_PATH` is a path that points to the location of your pseudopotential files. These must have the naming convention: `[element name].pot` to be detected. Currently SPARC only uses the psp8 format. Alternatively, rather than setting this environment variable you may pass the information in with the `pseudo_dir` argument.

~~~
#get sparc calculator
from sparc_core import SPARC
calc = SPARC(h=0.2) # a grid spacing grid must be entered.

#make atoms
from ase.build import bulk
atoms = bulk('Si',cubic=True) #cell must be rectangular, This feature is under development
atoms.set_calculator(calc)
atoms.get_potential_energy()
~~~

