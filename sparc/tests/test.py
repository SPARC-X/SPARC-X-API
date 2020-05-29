import os
from copy import deepcopy as dc
from sparc.sparc_core import SPARC
from sparc.ion import write_ion, read_ion
from ase.build import bulk, molecule
from ase.calculators.calculator import compare_atoms
from ase.visualize import view
from ase.constraints import FixAtoms, FixedLine, FixedPlane
import numpy as np

atoms = bulk('NaCl', crystalstructure = 'rocksalt', a = 5) * (2,1,1)
atoms.set_initial_magnetic_moments([1.1] * len(atoms))
l = SPARC(atoms = atoms, h = 0.1, label = 'in1', xc = 'GGA', spin_typ = 2)
l.write_input(atoms = atoms, h =0.1, spin_typ = 1)
print(l.parameters)
print(l.estimate_memory())

print('input writing functional')
# check reading and writing .ion files
atoms.set_constraint([FixAtoms([0]), FixedLine(1,[0,1,0]), FixedPlane(2,[1,0,0])])
write_ion(open('in1.ion','w'), atoms)
recovered_atoms = read_ion(open('in1.ion','r'))
assert compare_atoms(atoms, recovered_atoms) == []
print('.ion reading test passed')
#view(recovered_atoms)
#view(atoms)


calc = SPARC.read('sprc-calc')

print('reading test passed')

d = calc.todict()
print('to dict test passed')

