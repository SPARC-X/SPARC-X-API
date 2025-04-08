# Using SPARC
import numpy as np
from ase.build import molecule

from sparc.calculator import SPARC

water = molecule("H2O", vacuum=7)
water.pbc = [False, False, False]

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1, 1, 1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "ELEC_TEMP_TYPE": "fermi-dirac",
    "ELEC_TEMP": 116,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "RELAX_FLAG": 1,
}


def main():
    """Geometric optimization of a water molecule using
    ab initio SPARC in FileIO mode
    """
    water.calc = SPARC(**calc_params)
    energy = water.get_potential_energy()
    forces = water.get_forces()

    print("Energy:", energy)
    print("Max Force:", np.max(forces))


if __name__ == "__main__":
    main()
