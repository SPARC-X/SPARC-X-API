# Simple DFT workflows with SPARC-X-API

As documentation for [basic usage](../basic_usage.md) shows, replacing an
existing DFT-based workflow using SPARC-X-API calculator interface can
be as easy as simply swapping the calculator instance to `sparc.SPARC`
from other codes like VASP, QE or GPAW.

Here we show a simple example that calculates the equation of state
(EOS) of bulk aluminum and determine its optimal lattice constant,
adapted from GPAW's
[tutorial](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/lattice_constants/lattice_constants.html). The full python script for this example can be downloaded [here](https://raw.githubusercontent.com/SPARC-X/SPARC-X-API/refs/heads/master/examples/simple_examples/ex0-eos.py).

In
the GPAW tutorial, a GPAW calculator in Planewave (PW) mode is created
like follows:

```{code} python
from gpaw import GPAW, PW
calc = GPAW(mode=PW(ecut),
            xc='PBE',
            kpts=(8, 8, 8),
            basis='dzp',
            txt=f'Al-{ecut}.txt')
```

We can create a SPARC-X-API calculator using similar parameters. Note
that in real-space DFT, the parameter mesh spacing (`h`) controls the
convergence. To avoid large "egg-box effect" due to large mesh
spacing, we recommend to use a smaller `h` value. For demonstration
purpose a rather rough mesh spacing `h=0.25` (in Angstrom) and a 3x3x3
k-points are used.

```{code} python
import numpy as np
from ase.build import bulk
from ase.eos import calculate_eos

from sparc import SPARC


def main():
    # Al in conventional cell
    atoms = bulk("Al", cubic=True)
    calc = SPARC(h=0.25, kpts=(3, 3, 3), xc="pbe", directory="ex0-eos")
    vol = atoms.get_volume()
    atoms.calc = calc
    eos = calculate_eos(atoms, npoints=5, eps=0.05, trajectory="al-eos-sparc.traj")
    print("Original volume: Ang^3", vol)
    v, e, B = eos.fit()
    print("Fitted volume (Ang^3), energy (eV), modulus (eV/Ang^3)")
    print(v, e, B)
    a0 = v ** (1 / 3)
    print(f"Optimal cell length (cubic): {a0} Ang")
    atoms.set_cell([a0, a0, a0], scale_atoms=True)
    e0_sparc = atoms.get_potential_energy()
    print(f"Energy calculated by SPARC: {e0_sparc} eV")
    print(f"Energy diff {abs(e0_sparc - e)} eV")
    return


if __name__ == "__main__":
    main()
```

The output from the above example may look like:
```{raw}
Fitted volume (Ang^3), energy (eV), modulus (eV/Ang^3)
65.97840834969949 -253.07755156337953 2.9095110471623173
Optimal cell length (cubic): 4.040799280428726 Ang
Energy calculated by SPARC: -253.0552324051582 eV
Energy diff 0.02231915822133601 eV
```

```{note}
This example uses file I/O mode for demonstration purpose only. Consider choosing the [socket mode](../advanced_socket.md) if you need more flexibility.
```
