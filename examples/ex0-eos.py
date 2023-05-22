"""A simple example using sparc-dft-api to calculate the equation of state of bulk Al.

Example taken from GPAW's tutorial:
https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/lattice_constants/lattice_constants.html
"""
import numpy as np
from ase.build import bulk
from sparc import SPARC
from ase.eos import calculate_eos

def main():
    # Al in conventional cell
    atoms = bulk('Al', cubic=True)
    calc = SPARC(h=0.30, kpts=(3, 3, 3), xc="pbe", directory="ex0-eos")
    vol = atoms.get_volume()
    atoms.calc = calc
    eos = calculate_eos(atoms, npoints=5, eps=0.05, trajectory="al-eos-sparc.traj")
    print("Original volume: ", vol)
    v, e, B = eos.fit()
    print(v, e, B)
    return

if __name__ == "__main__":
	main()

