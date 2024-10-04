# Using SPARC
from sparc.calculator import SPARC
from ase.optimize import BFGS
from ase.build import molecule

water = molecule('H2O', vacuum=7)
water.pbc = [False,False,False]

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
    "SPIN_TYP": 0,
}

with SPARC(use_socket=True, **calc_params) as calc:
    water.calc = calc
    dyn = BFGS(water, trajectory = 'water-bfgs-opt.traj', logfile= 'opt.log')
    dyn.run(fmax=0.05)