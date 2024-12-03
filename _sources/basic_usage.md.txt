# Basic Usage

## Read / write SPARC files

In contrast to many other DFT codes, where the ASE I/O formats refer
to a single file, `SPARC-X-API` operates on the whole calculation
directory, also known as a "SPARC bundle". This API integrates
seamlessly with ASE, allowing for the automatic detection of the SPARC
file format:

- Reading from a SPARC bundle

```python
import sparc
from ase.io import read, write
atoms = read("test.sparc", index=-1, format="sparc")
```
*Note*: To read multiple output files from the same directory, e.g., SPARC.aimd, SPARC.aimd\_01, pass the keyword argument `include_all_files=True` to `read()`

- Writing a minimal SPARC bundle from atoms

```python
import sparc
from ase.io import read, write
from ase.build import Bulk
atoms = Bulk("Al") * [4, 4, 4]
atoms.write("test.sparc", format="sparc")
```

```{note}
You need to specify `format="sparc"` when using the `read` and `write` functions from `ase.io`, as automatic file extension detection doesn't work for directories.
```

For a deeper dive into the bundle I/O format, see [Advanced Topics](advanced_topics.md).

### JSON Schema for SPARC calculator

A recurring challenge of Python interfaces to DFT codes it the
inconsistencies between low-level codes (Fortran/C/C++) and outdated
upper-level APIs regarding parameter sets and default values. To
address this issue, `SPARC-X-API` handles DFT parameters through a
JSON schema translated from SPARC's LaTeX documentation.  Each release
of `SPARC-X-API` is linked with a specific version of the SPARC
source code, ensuring compatibility and consistency with the default
parameter set. The main driver of this feature is the
`sparc.api.SparcAPI` class.

If you've obtained the full SPARC [source
code](https://github.com/SPARC-X/SPARC), you can generate a copy of
the schema by the following code:
```bash
python -m sparc.docparser <sparc-source-code-root>/doc/.LaTeX
```
which produces a `parameters.json` file.

To learn more about the JSON schema design, please refer to [Advanced
Topics](advanced_topics.md).


### Calculator interface (File-IO mode)

`SPARC-X-API` offers a calculator interface based on file I/O that aligns with many
other ASE calculators.  If you've worked with ASE modules like `Vasp`,
`QuantumEspresso`, or `GPAW`, you'll find this package intuitive,
as shown in the following examples:

1. Single point calculation
```python
from sparc.calculator import SPARC
from ase.build import molecule
atoms = molecule("H2", cell=(10, 10, 10), pbc=True, directory="run_sp")
atoms.calc = SPARC(h=0.25)
atoms.get_potential_energy()
atoms.get_forces()
```

This example sets up a calculation for H2 atoms in a 10 Å x 10 Å x 10
Å PBC cell with default parameters (PBE exchange correlation
functional, and a grid spacing (`h`) of 0.25 Å). Note by calling
`atoms.get_forces`, the calculator will automatically sets the flags
for printing the forces.

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

This example sets up a calculation for a rattled Aluminum primitive
unit cell, calculate with PBE functional, grid spacing of 0.25 Å, and
3 x 3 x 3 k-point grid.  Optimization of ionic positions is handled
with SPARC's default LBFGS routine.

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

The power of `SPARC-X-API` is to combine single point `SPARC` calculations with advanced ASE optimizers, such as `BFGS`, `FIRE` or `GPMin`. Example 2 can be re-written as:

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

### Command-line tools

A simple command wrapper `sparc-ase` is provided to add
support of SPARC file formats to the `ase` cli tools. Simple
replace `ase [subcommand] [args]` with `sparc-ase [subcommand] [args]`
to access your SPARC bundle files as you would use for other file
formats.  As an example, use `sparc-ase gui path/to/your/bundle.sparc`
for the visualization of atomistic structures. Depending on the
bundle's contents, this could display individual atoms or multiple
images.

[Fig. 1](#fig-screenshot-sparc-ase) is a screenshot showing the usage of `sparc-ase gui` to visualize a
short MD trajectory.

(fig-screenshot-sparc-ase)=
```{figure} https://github.com/alchem0x2A/SPARC-X-API/assets/6829706/e72329ff-7194-4819-94f8-486ef2218844

Fig 1. A screenshot of the `sparc-ase` program
```

### Parameters and units used in `SPARC-X-API`

In the SPARC DFT code, all input parameters conventionally employ atomic units, such as Hartree and Bohr. Conversely, ASE objects (like `Atoms.positions`, `Atoms.cell`, `Atoms.get_potential_energy()`) utilize eV/Angstrom units.

When you set up a calculator as below:
```python
atoms.calc = SPARC(h=0.25, REFERENCE_CUTOFF=0.5, EXX_RANGE_PBE=0.16, **params)
```
inputs following ASE's convention (e.g., `h`) adopt eV/Angstrom units (thus the same setting can be applied to other DFT calculators),
On the other hand, all SPARC-specific parameters, which can often be recognized by their capitalized format (like `REFERENCE_CUTOFF`, `EXX_RANGE_PBE`), retain their original values consistent with their representation in the `.inpt` files.
The reasoning and details about unit conversion can be found in the [Rules for Input Parameters](#rule-input-param)  in Advanced Topics.


In order for `SPARC-X-API` to be compatible with other ASE-based DFT calculators,
there is a list of special parameters consistent with the ASE convention and uses Å / eV / GPa / fs
unit system:

| parameter name | meaning                         | example        | equivalent `SPARC` input |
|----------------|---------------------------------|----------------|--------------------------|
| `xc`           | Exchange-correlation functional | `xc=pbe` | `EXCHANGE_CORRELATION: GGA_PBE` |
|                |                                 | `xc=lda` | `EXCHANGE_CORRELATION: LDA_PZ` |
|                |                                 | `xc=rpbe` | `EXCHANGE_CORRELATION: GGA_RPBE` |
|                |                                 | `xc=pbesol` | `EXCHANGE_CORRELATION: GGA_PBEsol` |
|                |                                 | `xc=pbe0` | `EXCHANGE_CORRELATION: PBE0` |
|                |                                 | `xc=hf` | `EXCHANGE_CORRELATION: HF` |
|                |                                 | `xc=hse` or `xc=hse03` | `EXCHANGE_CORRELATION: HSE` |
|                |                                 | `xc=vdwdf1` or `xc=vdw-df` | `EXCHANGE_CORRELATION: vdWDF1` |
|                |                                 | `xc=vdwdf2` or `xc=vdw-df2` | `EXCHANGE_CORRELATION: vdWDF2` |
|                |                                 | `xc=scan` | `EXCHANGE_CORRELATION: SCAN` |
| `h`            | Real grid spacing    (Å)           | `h=0.2`        | `MESH_GRID: 0.38`  (in Bohr)         |
| `gpts`         | Explicit grid points |   `gpts=[10, 10, 10]` | `FD_GRID: 10 10 10` |
| `kpts`         | Kpoint mesh          |   `kpts=[3, 3, 3]`    | `KPOINT_GRID: 3 3 3` |
| `convergence`  | Dict of convergence criteria (see below) |  | |
|                | `energy`  eV/atom         | `convergence={"energy": 1e-4}` | `TOL_SCF: 3e-6` |
|                | `relax` (forces)  eV/Å            | `convergence={"relax": 1e-2}` | `TOL_RELAX: 2e-4` |
|                | `density` e/atom          | `convergence={`density`: 1e-6}`| `TOL_PSEUDOCHARGE: 1e-6` |

Users from other DFT codes can easily port their ASE codes to `SPARC-X-API` using the special parameters with minimal modification:

Example 1: VASP vs SPARC

```python
# Using VASP
from ase.calculators.vasp import Vasp
calc = Vasp(xc="rpbe", kpts=(9, 9, 9), directory="vasp-calc")
```
vs
```python
# Using SPARC
from sparc.calculator import SPARC
calc = SPARC(xc="rpbe", kpts=(9, 9, 9), directory="sparc-calc.sparc")
```

Example 2: GPAW (another real-space DFT code) vs SPARC
```python
# Using GPAW
from gpaw import GPAW
calc = GPAW(xc="PBE", kpts=(9, 9, 9), h=0.25, directory="gpaw-calc", convergence={"energy": 1.e-4})
```
vs
```python
# Using SPARC
from sparc.calculator import SPARC
calc = SPARC(xc="PBE", kpts=(9, 9, 9), h=0.25, directory="sparc-calc.sparc", convergence={"energy": 1.e-4})
```
