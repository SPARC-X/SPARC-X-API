from ase import Atoms
from ase.constraints import FixedPlane
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.units import fs

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
cons = [FixedPlane(i, [1, 0, 0]) for i in range(5)]
Ag_cluster.set_constraint(cons)

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
}

taut = 50 * fs
timestep = 2 * fs


def main():
    """MD simulation with implicit socket interface"""
    with SPARC(use_socket=True, **calc_params) as calc:
        Ag_cluster.calc = calc
        # dyn = NVTBerendsen(Ag_cluster, temperature_K=10, taut=taut, timestep=timestep, trajectory='Ag-cluster-md.traj')
        dyn = VelocityVerlet(
            Ag_cluster,
            dt=2.0 * fs,
            trajectory="Ag-cluster-verlet.traj",
            logfile="verlet-md.log",
        )
        dyn.run(10)


if __name__ == "__main__":
    main()
