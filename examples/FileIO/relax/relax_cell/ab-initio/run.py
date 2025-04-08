# Using SPARC
from ase.build import bulk

from sparc.calculator import SPARC

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [4, 4, 4],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "CALC_STRESS": 1,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "RELAX_FLAG": 2,
}


def main():
    atoms = bulk("Rh", crystalstructure="fcc", a=3.81)
    atoms.calc = SPARC(**calc_params)
    # Trigger the calculation by calling a property that will be available after the calculation
    stress = atoms.get_stress()
    # Get the final energy and volume from a single point calculation
    del calc_params["RELAX_FLAG"]
    if type(atoms) == list:
        atoms = atoms[-1]
    atoms.calc = SPARC(**calc_params)
    energy = atoms.get_potential_energy()
    volume = atoms.get_volume()

    print(
        "***********************************************************************************************"
    )
    print("v0 = {0} A^3\nE0 = {1} eV".format(volume, energy))


if __name__ == "__main__":
    main()
