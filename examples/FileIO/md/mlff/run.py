import numpy as np
from ase import Atoms

from sparc.calculator import SPARC

Ag_cluster = Atoms(
    "Ag5",
    positions=[
        (0.0, 2.6579, 0.9366),
        (0.0, -1.3587, -1.4045),
        (0.0, 0.0, 0.9358),
        (0.0, -2.6579, 0.9366),
        (0.0, 1.3587, -1.4045),
    ],
    pbc=(0, 0, 0),
)
Ag_cluster.set_cell([20.0, 24.0, 24.0])
Ag_cluster.center()

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [1, 1, 1],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
    "MD_FLAG": 1,
    "MD_METHOD": "NVK_G",
    "ION_TEMP": 10,
    "MD_NSTEP": 10,
    "MD_TIMESTEP": 2,
    "MLFF_FLAG": 1,
    "MLFF_INITIAL_STEPS_TRAIN": 3,
}


def main():
    """Running a simple MD calculation using SPARC-MLFF in FileIO mode"""
    Ag_cluster.calc = SPARC(label="Ag_cluster", **calc_params)
    energy = Ag_cluster.get_potential_energy()
    print("Energy: ", energy)
    forces = Ag_cluster.get_forces()
    print("Max Force: ", np.max(forces))


if __name__ == "__main__":
    main()
