# Using SPARC
from ase import Atoms, units
from ase.calculators.plumed import Plumed
from ase.md.nvtberendsen import NVTBerendsen

from sparc.calculator import SPARC

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

timestep = (
    2 * units.fs
)  # the units module contains the conversion factor to go from units.<unit> to ASE units via multiplication
ps = 1000 * units.fs

# Units conversion: everything will be in ASE units, so a time step is in ASE unit and needs to be interpreted as some multiple of sensical units
# Likewise, energy is in eV which is ~98.6 kJ/mol, so a single unit of energy is 98.6 kJ/mol which are the internal units of PLUMED
setup = [
    f"UNITS LENGTH=A TIME={1/ps} ENERGY={units.mol/units.kJ}",  # Set units to match desired properties
    # Calculate the center of mass of atoms 1-5
    "com: COM ATOMS=1-5",
    # Define the coordination number (C)
    "c: COORDINATION GROUPA=1-5 SWITCH={RATIONAL R_0=3.0 NN=8 MM=16} NOPBC",
    # Define the radius of gyration (R)
    "r: GYRATION TYPE=RADIUS ATOMS=1-5 NOPBC",
    # Compute CV1 and CV2 as orthogonal linear combinations of C and R
    "cv1: COMBINE ARG=c,r COEFFICIENTS=0.99715,-0.07534 PERIODIC=NO",
    "cv2: COMBINE ARG=c,r COEFFICIENTS=0.07534,0.99715 PERIODIC=NO",
    # Apply lower wall on CV1 at 5.0 with harmonic constant 10 eV
    f"LOWER_WALLS ARG=cv1 AT=5.0 KAPPA=10.0",
    # Apply lower wall on CV2 at 3.0 with harmonic constant 50 eV
    f"UPPER_WALLS ARG=cv2 AT=3.0 KAPPA=50.0",
    # Perform well-tempered metadynamics on CV1 and CV2
    f"METAD ARG=cv1,cv2 HEIGHT=0.3 PACE=500 SIGMA=0.3,0.03 GRID_MIN=0.0,0.0 GRID_MAX=10.0,5.0 GRID_BIN=500,500 BIASFACTOR=100 FILE=HILLS",
    # Print out the collective variables for monitoring
    "PRINT ARG=cv1,cv2,c,r FILE=COLVAR STRIDE=1",
]

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
# cons = [FixedPlane(i, [1, 0, 0]) for i in range(5)]
# Ag_cluster.set_constraint(cons)


def main():
    """Sampling via Plumed interface with SPARC's socket"""
    with SPARC(use_socket=True, **calc_params) as calc:
        Ag_cluster.calc = Plumed(
            calc=calc, input=setup, timestep=timestep, atoms=Ag_cluster, kT=0.00861733
        )  # 10 K in eV thermal energy units
        dyn = NVTBerendsen(
            Ag_cluster,
            timestep,
            temperature_K=0.00861733 / units.kB,
            taut=50 * units.fs,
            fixcm=False,
            trajectory="Ag-cluster-metadynamics.traj",
        )
        dyn.run(10)

        """
    # Restrain atoms 1-5 within 2.0 Å of the center of mass using upper walls
    "d1: DISTANCE ATOMS=1,com",
    "UPPER_WALLS ARG=d1 AT=2.0 KAPPA=100.0",
    "d2: DISTANCE ATOMS=2,com",
    "UPPER_WALLS ARG=d2 AT=2.0 KAPPA=100.0",
    "d3: DISTANCE ATOMS=3,com",
    "UPPER_WALLS ARG=d3 AT=2.0 KAPPA=100.0",
    "d4: DISTANCE ATOMS=4,com",
    "UPPER_WALLS ARG=d4 AT=2.0 KAPPA=100.0",
    "d5: DISTANCE ATOMS=5,com",
    "UPPER_WALLS ARG=d5 AT=2.0 KAPPA=100.0",
　　　"""


if __name__ == "__main__":
    main()
