"""A simple example using SPARC-X-API to calculate the equation of state of bulk Al.

Example taken from GPAW's tutorial:
https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/lattice_constants/lattice_constants.html
"""
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
