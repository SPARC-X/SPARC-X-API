# Using SPARC
from sparc.calculator import SPARC
from ase.eos import EquationOfState
from ase.build import bulk

calc_params = {
    "EXCHANGE_CORRELATION": "GGA_PBE",
    "KPOINT_GRID": [4,4,4],
    "MESH_SPACING": 0.35,
    "TOL_SCF": 0.0001,
    "MAXIT_SCF": 100,
    "CALC_STRESS": 1,
    "PRINT_RESTART_FQ": 10,
    "PRINT_ATOMS": 1,
    "PRINT_FORCES": 1,
    "SPIN_TYP": 0,
}

LC = [3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0]
energies = []
volumes = []
with SPARC(use_socket=True, **calc_params) as calc:
    for a in LC:
        atoms = bulk('Rh', crystalstructure='fcc', a = a)
        atoms.calc = calc
        volumes.append(atoms.get_volume())
        energies.append(atoms.get_potential_energy())

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()

print('v0 = {0} A^3\nE0 = {1} eV\nB  = {2} eV/A^3'.format(v0, e0, B))

eos.plot('Rh-fcc-eos.png')