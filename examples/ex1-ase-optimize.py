"""A simple example using SPARC to optimize a NH3 molecule using:
1) SPARC internal LBFGS routine
2) SPARC single point + ASE LBFGS


"""
from sparc import SPARC

import numpy as np
from ase.build import molecule
from ase.optimize.lbfgs import LBFGS
from ase.constraints import FixAtoms

nh3 = molecule("NH3", cell=(6, 6, 6))
# Fix the N center
nh3.constraints = [FixAtoms([0])]
nh3.rattle()


def optimize_sparc_internal():
    atoms = nh3.copy()
    calc = SPARC(h=0.25, kpts=(1, 1, 1), xc="pbe",
                 convergence={"forces": 0.02},
                 relax_flag=True, print_relaxout=True,
                 relax_method="LBFGS",
                 directory="ex1-sparc")
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
    calc = SPARC(h=0.25, kpts=(1, 1, 1), xc="pbe", print_force=True,
                 directory="ex1-ase")
    atoms.calc = calc
    opt = LBFGS(atoms)
    opt.run(fmax=0.02)
    # Number of ionic steps in case calc.get_number_of_ionic_steps not implemented
    # nsteps = len(calc.raw_results["geopt"])
    # print("SPARC internal LBFGS:")
    # print(f"Final energy: {e_fin} eV")
    # print(f"Final fmax: {np.max(np.abs(f_fin))} eV/Ang")
    # print(f"N steps: {nsteps}")
    
    
    # # Al in conventional cell
    # atoms = bulk('Al', cubic=True)
    # calc = SPARC(h=0.25, kpts=(3, 3, 3), xc="pbe", directory="ex0-eos")
    # vol = atoms.get_volume()
    # atoms.calc = calc
    # eos = calculate_eos(atoms, npoints=5, eps=0.05, trajectory="al-eos-sparc.traj")
    # print("Original volume: Ang^3", vol)
    # v, e, B = eos.fit()
    # print("Fitted volume (Ang^3), energy (eV), modulus (eV/Ang^3)")
    # print(v, e, B)
    # a0 = v ** (1 / 3)
#     # print(f"Optimal cell length (cubic): {a0} Ang")
#     # atoms.set_cell([a0, a0, a0], scale_atoms=True)
#     # e0_sparc = atoms.get_potential_energy()
#     # print(f"Energy calculated by SPARC: {e0_sparc} eV")
#     # print(f"Energy diff {abs(e0_sparc - e)} eV")
#     return

if __name__ == "__main__":
    optimize_sparc_internal()

