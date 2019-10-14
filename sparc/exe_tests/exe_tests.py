from sparc.sparc import SPARC
from ase.build import molecule, bulk
import os

molecule_test = molecule('H2O')
molecule_test.set_cell([6,6,6])
molecule_test.center()

#### basic_funcitonality
# Direlecht boundary conditions, SCF_TOL inputs, medium grids
# Absolute coordinates

calc = SPARC(atoms = molecule_test)
calc.write_input(atoms = molecule_test, EXCHANGE_CORRELATION = 'LDA_PZ', 
                 label = 'H2O', scaled = False,
                 h = 0.2, TOL_SCF = 1e-4,
                 BOUNDARY_CONDITION = 1, KPOINT_GRID = (1,1,1))

os.system('mv H* O* mol_dirlchet')


# fermi-dirac smearing, periodic boundary conditions
# coarse grids, K-points, density mixing, kerker preconditioner, d block
bulk_test = bulk('Cu', cubic = True)

calc = SPARC(atoms = molecule_test,BOUNDARY_CONDITION = 2, KPOINT_GRID = (2,2,2))
calc.write_input(atoms = bulk_test, label = 'Cu', SMEARING = 0.01,
                 ELEC_TEMP_TYPE = 'fermi-dirac',
                 MIXING_VARIABLE = 'density', MIXING_PRECONDITIONER = 'kerker',
                 h = 0.3, TOL_SCF = 1e-4,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (2,2,2))

os.system('mv Cu* bulk_pbc_kpts')


#### Relaxation

# Relaxation, LBFGS, fine grids, minimal accuracy, GGA

opt_test = molecule('H2')
opt_test.set_cell([5,5,5])
opt_test.center()

calc = SPARC(atoms = opt_test)
calc.write_input(atoms = opt_test, label = 'H2',
                 EXCHANGE_CORRELATION = 'GGA_PBE',
                 h = 0.1, TOL_SCF = 1e-3,
                 RELAX_FLAG = 1,
                 #ACCURACY = 'minimal'
                 RELAX_METHOD = 'LBFGS',
                 L_HISTORY = '18',
                 L_FINIT_STP = 0.004,
                 L_MAXMOV = 0.21,
                 BOUNDARY_CONDITION = 1, KPOINT_GRID = (1,1,1))

os.system('mv H* LBFGS_relax_test')

opt_test = molecule('O2')
opt_test.set_cell([5,5,5])
opt_test.center()

# scf energy accuracy check, NLCG check, electronic temperature 
# kerker with potential mixing, LDA pw, 2D periodic

calc = SPARC(atoms = opt_test)
calc.write_input(atoms = opt_test, label = 'O2',
                 #ACCURACY = 'medium',
                 EXCHANGE_CORRELATION = 'LDA_PW',
                 ELEC_TEMP = 1100,
                 SCF_ENERGY_ACC = 0.002,
                 MIXING_VARIABLE = 'potential',
                 MIXING_PRECONDITIONER = 'kerker',
                 RELAX_FLAG = 1,
                 RELAX_METHOD = 'NLCG',
                 BOUNDARY_CONDITION = 3, KPOINT_GRID = (1,1,1))


os.system('mv O* NLCG_relax_test')

# scf force accuracy check, FIRE check, mixing parameter check
# relax printing, relax tol, 1D periodic 
opt_test = molecule('P2')
opt_test.set_cell([5,5,5])
opt_test.center()

calc = SPARC(atoms = opt_test)
calc.write_input(atoms = opt_test, label = 'P2',
                 #ACCURACY = 'medium',
                 SCF_FORCE_ACC = 0.0002,
                 MIXING_PARAMETER = 0.35,
                 TOL_RELAX = 2e-3,
                 RELAX_FLAG = 1,
                 PRINT_RELAXOUT = 0,
                 RELAX_METHOD = 'FIRE',
                 BOUNDARY_CONDITION = 4, KPOINT_GRID = (1,1,1))


os.system('mv P* FIRE_relax_test')

#### Stress/Pressure in Bulks


# Stress, Pressure, density with no preconditioner
# high accuracy

stress_test = bulk('Li', cubic = True)

calc = SPARC(atoms = molecule_test)
calc.write_input(atoms = stress_test, label = 'Li',
                 MIXING_VARIABLE = 'density', MIXING_PRECONDITIONER = 'none',
                 CALC_STRESS = 1, CALC_PRES = 1,
                 h = 0.3, TOL_SCF = 1e-6,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (4,4,4))


os.system('mv Li* stress_pressure_test')


#### print functions

print_test = bulk('La')

# all print functions, f block element, non-orthogonal unit cells
calc = SPARC(atoms = molecule_test)
calc.write_input(atoms = print_test, label = 'La',
                 MIXING_VARIABLE = 'density', MIXING_PRECONDITIONER = 'kerker',
                 CALC_STRESS = 1, CALC_PRES = 1,PRINT_FORCES = 1, PRINT_ATOMS = 1,
                 PRINT_EIGEN = 1, PRINT_RESTART_FQ = 1, PRINT_RESTART =1,
                 PRINT_VELS = 1,
                 h = 0.35, TOL_SCF = 1e-3,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (2,2,2))


os.system('mv La* print_test')


#### MD

print_test = bulk('C', cubic = True)


# MD, NVT, FD order

calc = SPARC(atoms = molecule_test)
calc.write_input(atoms = print_test, label = 'C',
                 h = 0.35, TOL_SCF = 1e-3,
                 MD_FLAG = 1,MD_METHOD = 'NVT_NH',
                 FD_ORDER = 10,
                 ION_TEMP = 100,
                 MD_NSTEP = 5,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (1,1,1))

os.system('mv C* MD_NVT_test')

# MD, NVE, Q mass, Beta 

print_test = bulk('Na', cubic = True)

calc = SPARC(atoms = molecule_test)
calc.write_input(atoms = print_test, label = 'Na',
                 QMASS = 1.1,
                 BETA = 1100,
                 h = 0.35, TOL_SCF = 1e-3,
                 MD_FLAG = 1,MD_METHOD = 'NVE',
                 ION_TEMP = 100,
                 MD_NSTEP = 5,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (1,1,1))

os.system('mv Na* MD_NVE_test')

### Tols

# strange tol combinations

bulk_test = bulk('Mg', cubic = True)

calc = SPARC(atoms = molecule_test,BOUNDARY_CONDITION = 2, KPOINT_GRID = (2,2,2))
calc.write_input(atoms = bulk_test, label = 'Mg',
                 MIXING_VARIABLE = 'density', MIXING_PRECONDITIONER = 'kerker',
                 TOL_POISSON = 1e-4, TOL_LANCZOS = 5e-3,
                 h = 0.15, TOL_SCF = 1e-4,
                 BOUNDARY_CONDITION = 2, KPOINT_GRID = (2,2,2))

os.system('mv Mg* tol_test')


