import os
import numpy as np
import warnings

from utilities import h2gpts
from ase.units import Bohr, Hartree
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import CalculatorError, CalculatorSetupError
from ase.calculators.calculator import EnvironmentError, InputError
from ase.calculators.calculator import CalculationFailed, SCFError, ReadError
from ase.calculators.calculator import PropertyNotImplementedError, PropertyNotPresent

from ion import write_ion

required_inputs = ['PSEUDOPOTENTIAL_FILE',
                    'CELL','EXCHANGE_CORRELATION',
                    'FD_GRID','PSEUDOPOTENTIAL_LOCAL',
                    'KPOINT_GRID', 'LATVEC']

default_parameters = {
            # 'label': 'sprc-calc',
            # 'calculation_directory':'sprk.log',

            'BOUNDARY_CONDITION': 2,
            'LATVEC': None,
            'EXCHANGE_CORRELATION': 'LDA_PZ',  # 'LDA'
            'KPOINT_GRID': (1, 1, 1),
            'MIXING_PARAMETER': 0.30,
            'CHEN_DEGREE': 20,
            'NSTATES': None,
            'SMEARING': None,
            'MAXIT_SCF': 100,
            'BETA': 1000,
            'ELEC_TEMP': None,
            'CALC_STRESS': None,
            'CALC_PRES': None,
            'TOL_SCF': 1.00E-05,
            'TOL_POISSON': 1.00E-06,
            'TOL_LANCZOS': 1.00E-02,
            'TOL_PSEUDOCHARGE': 1.00E-08,
            'TWTIME': 999999999.000000,
            'MIXING_PARAMETER': 0.30,
            'MIXING_HISTORY': 7,
            'MIXING_VARIABLE': None,
            'MIXING_PRECOND': None,
            'PULAY_FREQUENCY': 1,
            'PULAY_RESTART': 0,
            'REFERENCE_CUTOFF': 0.50,
            'RHO_TRIGER': 3,
            'PRINT_FORCES': 0,
            'PRINT_ATOMS': 0,
            'PRINT_EIGEN': 0,
            'PRINT_DENSITY': 0,
            'PRINT_RESTART_FQ': 1,
            'PRINT_RESTART': 1,
            'PSEUDOPOTENTIAL_LOCAL': None,
            'PSEUDOPOTENTIAL_FILE': None,
            'OUTPUT_FILE': None,
            'CELL': None,
            'FD_GRID': None,
            'FD_ORDER': 12,
            'ELEC_TEMP': 315.775131,
            'ELEC_TEMP_TYPE': None,
            'CHEB_DEGREE': 25,
            #'NTYPES': None,
            'NP_KPOINT_PARAL': None,
            'NP_BAND_PARAL': None,
            'NP_DOMAIN_PARAL': None,
            'NP_DOMAIN_PHI_PARAL': None,

            'TOL_RELAX': 1.00E-03,
            'PRINT_RELAXOUT': 0,
            'RELAX_FLAG': 0,
            'RELAX_METHOD': None,
            'RELAX_MAXITER': 300,
            'RELAX_NITER': 300,
            'NLCG_sigma': 0.500000,
            'L_HISTORY': 20,
            'L_FINIT_STP': 0.005000,
            'L_MAXMOV': 0.200000,
            'L_AUTOSCALE': 1,
            'L_LINEOPT': 1,
            'L_ICURV': 1.000000,

            'MD_FLAG': None,
            'MD_METHOD': None,
            'MD_TIMESTEP': None,
            'MD_NSTEP': None,
            'PRINT_RESTART_FQ': None,
            'RESTART_FLAG': None,
            'ION_TEMP': None,
            'ION_ELEC_EQT': None,

                        }

equivalencies = {
            'xc': 'EXCHANGE_CORRELATION',
            'kpts': 'KPOINT_GRID',
            'nbands': 'NSTATES',
            'h': None,

            'gpts': 'FD_GRID'

            }


class SPARC(FileIOCalculator):


    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='sprc-calc', atoms=None, command=None, directory = '.',
                 **kwargs):
        FileIOCalculator.__init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, command=None, directory = directory, **kwargs)

        # setting up label
        ## TODO: figure out how to do this with ASE calculator default classes
        if self.directory != '.' and '/' in label:
            raise CalculatorSetupError('cannot set both directory and input `/`'
                                       ' in the label name')
        elif '/' in label:
            directory, label = label.split('/')
        self.directory = directory
        self.label = label

        # check that all arguements are legitimate arguements
        for arg in kwargs.keys():
            if arg.upper() not in default_parameters.keys() and arg not in equivalencies.keys():
                raise InputError('the argument {} was not found in the list of '
                                 'allowable arguments'.format(arg))

        FileIOCalculator.set(self, **kwargs)
        self.atoms = atoms

    def write_input(self, atoms = None, **kwargs):
        if atoms == None:
            if self.atoms == None:
                raise InputError('An atoms object must be provided or the '
                                 'calculator object must have atoms attached to'
                                 ' write an input file')
            atoms = self.atoms
        FileIOCalculator.write_input(self, atoms)

        f = open(os.path.join(self.directory, self.label + '.inpt'), 'w')

        # make all kwargs upper case
        for arg in kwargs:
            if arg.upper() in default_parameters:
                kwargs[arg.upper()] = kwargs.pop(arg)

        ############## Begin writing the file ####################
        f.write('# Input File Generated By SPARC ASE Calculator #\n')

        # deal with the finite differnce grid
        if 'h' in kwargs and 'FD_GRID' in kwargs:
            raise CalculatorSetupError('You cannot specify a grid spacing (h)'
                                       'and input the FD_GRID input arguement at'
                                       ' the same time')
        if 'h' not in kwargs and 'FD_GRID' not in kwargs:
            warnings.warn('neither a grid spacing (h) nor a finite difference '
                          'grid (FD_GRID) has been specified, this is not ideal.'
                          ' A default value of h = 0.15 angstrom has been inserted.')
            kwargs['h'] = 0.15
        
        if 'h' in kwargs:
           fd_grid = h2gpts(kwargs['h'], atoms.cell) 
        if 'FD_GRID' in kwargs:
            if type(kwargs['FD_GRID']) == str:
                fd_grid = kwargs['FD_GRID'].split()
                if len(fd_grid) != 3:
                    raise InputError('If a string is passed in for the FD_GRID flag it'
                                     ' must be able to split into a list of dimension 3'
                                     ' with .split()')
            fd_grid = list(kwargs['FD_GRID'])
            if len(fd_grid) != 3:
                raise InputError('if an iterable type is used for the FD_GRID flag, it'
                                 ' must be have dimension 3')
        f.write('{} {} {}\n'.format(*fd_grid))

        # Deal with the unit cell

        # this just saves a few levels of logic below
        if 'LATVEC' not in kwargs:
            kwargs['LATVEC'] = None
        if 'CELL' not in kwargs:
            kwargs['CELL'] = None

        # Scold/warn the user about using the CELL/LATVEC arguments
        if  kwargs['LATVEC'] is not None and  kwargs['CELL'] is None:
            raise InputError('If LATVEC is input, you also must provide the CELL argument')

        if  kwargs['LATVEC'] is None and  kwargs['CELL'] is not None:
            warnings.warn('The CELL argument was entered, but not the LATVEC argument.'
                          ' The unit cell in the input atoms object will be ignored and'
                          ' the cell will be assumed to be orthogonal')

        # Deal with CELL and LATVEC inputs
        if kwargs['CELL'] is not None:
            # If there's no LATVEC input, just assume it's an orthogonal unit cell
            if kwargs['LATVEC'] is None:
                    lattice = np.eye(3)
            elif type(kwargs['LATVEC']) == str:
                kwargs['LATVEC'] = kwargs['LATVEC'].split()
                kwargs['LATVEC'] = [float(a) for a in kwargs['LATVEC']]
                if len(kwargs['LATVEC']) != 9:
                    raise InputError('The value of LATVEC must have 9 elements (3x3)')
                lattice = np.array(kwargs['LATVEC'])
                lattice = np.reshape(lattice, (3,3))
            else:
                # Decipher the LATVEC input
                try:  # try to deal with numpy arrays by making them lists
                    kwargs['LATVEC'] = list(kwargs['LATVEC'])
                except:
                    pass
                # Check that LATVEC is a native iterable
                if type(kwargs['LATVEC']) != list:
                    raise InputError('LATVEC must be entered as a list, tuple, or '
                                     'space and linebreak separated string')
                # Deal with the simple case of lists or tuples
                else:
                    lattice = np.array(kwargs['LATVEC'])
                    lattice = np.reshape(lattice,(3,3))

            # Decipher CELL input
            if type(kwargs['CELL']) not in [str,list,tuple]:
                raise InputError('CELL must be entered as a list, tuple, or space'
                                 ' separated string')
            if type(kwargs['CELL']) == str:
                kwargs['CELL'] = kwargs['CELL'].split()
            kwargs['CELL'] = [float(a) for a in kwargs['CELL']]
            cell = np.array(kwargs['CELL'])
            if len(cell) != 3:
                raise InputError('The value of CELL must have 3 elements')
            atoms.cell = (lattice.T * cell).T * Bohr

        # If LATVEC and CELL aren't provided by the user, use the cell of the
        # atoms object to write the input (this is the prefered way)
        if round(float(np.linalg.norm(atoms.cell)), 1) == 0:
            f.write('CELL:')
            cell = np.eye(3) * (np.max(atoms.positions, axis=0) + (6, 6, 6))
            atoms.set_cell(cell)
            for cell_param in atoms.get_cell_lengths_and_angles()[0:3]:
                f.write(' ' + str((cell_param) / Bohr))
            atoms.center()
            f.write('\n')
        else:
            f.write('CELL:')
            for length in atoms.get_cell_lengths_and_angles()[:3]:
                f.write(' ' + str(length / Bohr))
            f.write('\nLATVEC:')
            for cell_vec in atoms.cell:
                lat_vec = cell_vec / np.linalg.norm(cell_vec)
                f.write('\n')
                for element in lat_vec:
                    f.write(str(element) + ' ')
            f.write('\n')

        # kpoints
        if 'kpts' in kwargs:
            if 'KPOINT_GRID' in kwargs:
                warnings.warn('both kpts and KPOINT_GRID were input, defaulting'
                              ' to kpts')
            kwargs['KPOINT_GRID'] = kwargs['kpts']
        if 'KPOINT_GRID' in kwargs:
            if kwargs['KPOINT_GRID'] is not None:
                f.write('KPOINT_GRID: ')
                if len(kwargs['KPOINT_GRID']) == 3:
                    for kpoint in kwargs['KPOINT_GRID']:
                        if type(kpoint) != int:
                            raise InputError('when KPOINT_GRID is entered as an'
                                             ' iterable, the values must be in the'
                                             ' integer type (i.e. (4,4,4))')
                        f.write(str(kpoint) + ' ')
                elif type(kwargs['KPOINT_GRID']) == str:
                    if len(kwargs['KPOINT_GRID'].split()) != 3:
                        raise InputError('when KPOINT_GRID is entered as a string, it'
                                        ' must have 3 elements separated by spaces '
                                        '(i.e. \'4 4 4\')')
                    f.write(kwargs['KPOINT_GRID'])
                elif len(kwargs['KPOINT_GRID']) != 3 or type(kwargs['KPOINT_GRID']) is not str:
                    raise InputError('KPOINT_GRID must be either a length 3 object'
                                     ' (i.e. (4,4,4)) or a string (i.e. \'4 4 4 \')')
                f.write('\n')

        # xc should be put in separately
        if 'xc' in kwargs:
            kwargs['EXCHANGE_CORRELATION'] = kwargs['xc']
            f.write('EXCHANGE_CORRELATION: ' +  
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
        else:
            if 'EXCHANGE_CORRELATION' not in kwargs:
                f.write('EXCHANGE_CORRELATION: ' +  
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
            else:
                f.write('EXCHANGE_CORRELATION: ' +  # Note the Miss-spelling
                        kwargs['EXCHANGE_CORRELATION'] + '\n')
        
        # turn on the flags for these print settings by default
        # it generally doesn't make sense to have these off
        for flag in ['PRINT_FORCES','PRINT_ATOMS']:
            if flag not in kwargs:
                f.write('{}: {}\n'.format(flag,'1'))

        # all non-required inputs
        for arg, value in kwargs.items():
            if arg not in required_inputs and arg not in equivalencies:
                f.write('{}: {}\n'.format(arg, value))
        f.close()

        # make the atomic inputs (.ion) file
        write_ion(open(self.label + '.ion','w'),
                  atoms, pseudo_dir = os.environ['SPARC_PSP_PATH'],
                  scaled = True)


    def setup_parallel_env(self):
        """
        just sets up the environment to have roughly optimal parallelization
        environment variables
        """
        if 'PBS_NP' in os.environ:
            if float(os.environ['PBS_NP']) > 1:
                os.environ['MV2_ENABLE_AFFINITY'] = '1'
                os.environ['MV2_CPU_BINDING_POLICY'] = 'bunch'
