# Using SPARC
from sparc.calculator import SPARC
from ase.build import molecule
import numpy as np

water = molecule('H2O', vacuum=7)
water.pbc = [False,False,False]
water.set_initial_magnetic_moments([0.1, 0.1, 0.1])

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1,1,1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "ELEC_TEMP_TYPE": "fermi-dirac",
    "ELEC_TEMP": 116,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 1,
}

with SPARC(use_socket=True, **calc_params) as calc:
    water.calc = calc
    print("Initial magnetic moments before calculation call:\n", water.get_initial_magnetic_moments())
    # energy = water.get_potential_energy()
    # forces = water.get_forces()
    net_magmom = water.get_magnetic_moment()
    magmoms = water.get_magnetic_moments()
# print('***********************************************************************************************')
# print("Energy: {}\nMax Force: {}".format(energy, np.max(abs(forces)) ) )
print('*'*100)
print("Net magnetic moment: ", net_magmom)
print("Atomic magnetic moments:\n", magmoms)