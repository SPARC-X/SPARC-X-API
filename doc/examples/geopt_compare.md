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
