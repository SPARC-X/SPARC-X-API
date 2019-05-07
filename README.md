# sparc-x_ase

sparc-x_ase is an ASE based python wrapper for the density functional theory (DFT) code SPARC. SPARC is under heavy software development, so much of this wrapper is likely to change. However, the wrapper is functional for single point energy calculations.

## Installation:
The following command should install the package:

~~~
pip install git+https://github.com/SPARC-X/pysparcx/
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

Currently installation is barebones, simply clone this repository (`git clone https://github.gatech.edu/bcomer3/sparc-x_ase.git`) and add the folder to your pythonpath (`export PYTHONPATH=[were you cloned the repository]:$PYTHONPATH`)


## Usage

To enable the SPARC calculator to work you need to tell it what command to use to run sparc and where your pseudopotentials are. These should be stored in the environment variables `$ASE_SPARC_COMMAND` and `$SPARC_PSP_PATH` respectively. e.g.:

`export ASE_SPARC_COMMAND=[command to run SPARC]`

`export SPARC_PSP_PATH=[location of pseudopotentials]`

The `ASE_SPARC_COMMAND` variable is somewhat tricky to set up. All the information about parallelization should be included where appropriate in the command (i.e. `mpirun -np X`.) SPARC takes in a flag `-name` followed by the prefix of the input file names. This wrapper deals with the prefixes completely internally, but the word `PREFIX` should be put in after the `-name` flag for SPARC. Here is an example command if the `sparc` executable is in your `$PATH` environment variable:

`mpirun -np $PBS_NP sparc -name PREFIX`

Once those environment varibles are in place, the SPARC ASE calculator works like any other ASE calculator. It must be imported, instantiated, and called here is some example code for calculating bulk Si:

~~~
#get sparc calculator
from sparc import SPARC
calc = SPARC()

#make atoms
from ase.build import bulk
atoms = bulk('Si',cubic=True) #cell must be rectangular, This feature is under development
atoms.set_calculator(calc)
atoms.get_potential_energy()
~~~

SPARC is still under development, so features are somewhat limited.
