"""A simple example using SPARC to optimize a NH3 molecule using:
1) SPARC internal LBFGS routine
2) SPARC single point + ASE LBFGS


"""
import numpy as np
from ase.build import molecule
from ase.constraints import FixAtoms

# from ase.optimize.lbfgs import LBFGS
from ase.optimize.bfgs import BFGS

from sparc import SPARC

nh3 = molecule("NH3", cell=(8, 8, 8), pbc=True)
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
    # breakpoint()
    e_fin = atoms.get_potential_energy()
    f_fin = atoms.get_forces()
    # Number of ionic steps in case calc.get_number_of_ionic_steps not implemented
    nsteps = len(calc.raw_results["geopt"])
    print("SPARC internal LBFGS:")
    print(f"Final energy: {e_fin} eV")
    print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    print(f"N steps: {nsteps}")


def optimize_ase_lbfgs():
    atoms = nh3.copy()
    calc = SPARC(
        h=0.25, kpts=(1, 1, 1), xc="pbe", print_forces=True, directory="ex1-ase"
    )
    atoms.calc = calc
    opt = BFGS(atoms)
    # breakpoint()
    opt.run(fmax=0.02)
    e_fin = atoms.get_potential_energy()
    f_fin = atoms.get_forces()
    nsteps = opt.nsteps
    print("ASE LBFGS")
    print(f"Final energy: {e_fin} eV")
    print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    print(f"N steps: {nsteps}")


if __name__ == "__main__":
    optimize_sparc_internal()
    optimize_ase_lbfgs()
