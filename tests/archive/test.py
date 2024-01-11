import os
from copy import deepcopy as dc

import numpy as np
from ase.build import bulk, molecule
from ase.calculators.calculator import CalculatorSetupError, compare_atoms
from ase.constraints import FixAtoms, FixedLine, FixedPlane
from ase.visualize import view

from sparc.ion import read_ion, write_ion
from sparc.sparc_core import SPARC

atoms = bulk("NaCl", crystalstructure="rocksalt", a=5) * (2, 1, 1)
atoms.set_initial_magnetic_moments([1.1] * len(atoms))
l = SPARC(atoms=atoms, h=0.1, label="in1", xc="GGA", spin_typ=1, pseudo_dir=".")
l.write_input(atoms=atoms, h=0.1, spin_typ=1, pseudo_dir=".", SPIN_TYP=1)

print("input writing functional")
# check reading and writing .ion files
atoms.set_constraint([FixAtoms([0]), FixedLine(1, [0, 1, 0]), FixedPlane(2, [1, 0, 0])])
write_ion(open("in1.ion", "w"), atoms, pseudo_dir=".")
recovered_atoms = read_ion(open("in1.ion", "r"))
assert compare_atoms(atoms, recovered_atoms) == []

calc = SPARC(atoms=atoms, pseudo_dir=".", SPIN_TYP=1)
try:  # check that `h` is required
    calc.write_input()
    raise Exception("test failed")
except CalculatorSetupError:
    pass

print(".ion reading test passed")

calc = SPARC.read("read_input")

print("reading test passed")

d = calc.todict()
print("to dict test passed")
