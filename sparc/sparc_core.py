import os
import numpy as np
import warnings
import inspect
import re
import subprocess
import json
from collections import OrderedDict

from ase.units import Bohr, Hartree, fs, GPa
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import CalculatorError, CalculatorSetupError
from ase.calculators.calculator import EnvironmentError, InputError
from ase.calculators.calculator import CalculationFailed, SCFError, ReadError
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.calculator import PropertyNotPresent
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import compare_atoms
from ase.io.trajectory import Trajectory
from ase.atoms import Atoms
from ase.atom import Atom
from ase.io.jsonio import encode

from .ion import write_ion, read_ion


special_inputs = ['PSEUDOPOTENTIAL_FILE',
                  'CELL', 'EXCHANGE_CORRELATION',
                  'FD_GRID', 'PSEUDOPOTENTIAL_LOCAL',
                  'pseudo_dir', 'KPOINT_GRID', 'LATVEC',
                  'copy_psp', 'BC']

default_parameters = {
    # 'label': 'sprc-calc',
    # 'calculation_directory':'sprk.log',

    'BOUNDARY_CONDITION': 2,
    'BC': 'P P P',
    'LATVEC': None,
    'EXCHANGE_CORRELATION': 'LDA_PZ',  # 'LDA'
    'MESH_SPACING': 0.4,
    'KPOINT_GRID': (1, 1, 1),
    'KPOINT_SHIFT': None,
    'MIXING_PARAMETER': 0.30,
    'CHEN_DEGREE': 20,
    'NSTATES': None,
    'SMEARING': None,
    'MAXIT_SCF': 100,
    'MINIT_SCF': 3,
    'MAXIT_POISSON': 3000,
    'POISSON_SOLVER': None,
    'BETA': 1000,
    'TOL_PRECOND': None,
    'PRECOND_KERKER_KTF': None,
    'ELEC_TEMP': None,
    'CALC_STRESS': None,
    'CALC_PRES': None,
    'TOL_SCF': 1.00E-05,
    'TOL_POISSON': 1.00E-06,
    'TOL_LANCZOS': 1.00E-02,
    'TOL_PSEUDOCHARGE': 1.00E-08,
    'TOL_RELAX_CELL': None,
    'SCF_ENERGY_ACC': None,
    'TWTIME': 999999999.000000,
    'MIXING_PARAMETER': 0.30,
    'MIXING_HISTORY': 7,
    'MIXING_VARIABLE': None,
    'MIXING_PRECOND': None,
    'PRECOND_KERKER_THRESH': None,
    'PULAY_FREQUENCY': 1,
    'PULAY_RESTART': 0,
    'REFERENCE_CUTOFF': 0.50,
    'RHO_TRIGGER': 3,
    'VERBOSITY': 1,
    'PRINT_FORCES': 0,
    'PRINT_ATOMS': 0,
    'PRINT_EIGEN': 0,
    'PRINT_DENSITY': 0,
    'PRINT_RESTART_FQ': 1,
    'PRINT_RESTART': 1,
    'PRINT_VELS': 1,
    'PSEUDOPOTENTIAL_LOCAL': None,
    'PSEUDOPOTENTIAL_FILE': None,
    'OUTPUT_FILE': None,
    'CELL': None,
    'FD_GRID': None,
    'FD_ORDER': 12,
    'ELEC_TEMP': 315.775131,
    'ELEC_TEMP_TYPE': None,
    'CHEB_DEGREE': 25,
    'CHEFSI_BOUND_FLAG': None,
    'FIX_RAND': None,

    'SPIN_TYP': 1,

    # 'NTYPES': None,
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
    'RELAX_MAXDILAT': None,
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
    'PRINT_MDOUT': None,
    'RESTART_FLAG': None,
    'ION_TEMP': None,
    'ION_ELEC_EQT': None,
    'ION_VEL_DSTR': None,
    'ION_TEMP_END': None,
    'QMASS': None,

}

equivalencies = {
    'xc': 'EXCHANGE_CORRELATION',
    'kpts': 'KPOINT_GRID',
            'nbands': 'NSTATES',
            'h': None,

            'gpts': 'FD_GRID'

}


class SPARC(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'stress']
    all_changes = ['positions', 'numbers', 'cell', 'pbc',
                   'initial_charges', 'initial_magmoms']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='sprc-calc', atoms=None, command=None, directory='.',
                 **kwargs):
        FileIOCalculator.__init__(self, restart=None, ignore_bad_restart_file=False,
                                  label=None, atoms=None, command=None, directory=directory, **kwargs)

        # setting up label
        if self.directory != '.' and '/' in label:
            raise CalculatorSetupError('cannot set both directory and input `/`'
                                       ' in the label name')
        elif '/' in label:
            directory, label = label.split('/')
        self.directory = directory
        self.label = label
        self.prefix = self.label
        self.results = {}

        FileIOCalculator.set(self, **kwargs)
        self.atoms = atoms

    def write_input(self, atoms=None, scaled=True,
                    **kwargs):

        if atoms is None:
            if self.atoms is None:
                raise InputError('An atoms object must be provided or the '
                                 'calculator object must have atoms attached to'
                                 ' write an input file')
            atoms = self.atoms
        if kwargs == {}:
            kwargs = self.parameters.copy()

        FileIOCalculator.write_input(self, atoms)

        # TODO: think about if this next conditional is a good idea
        if 'label' in kwargs:
            self.label = kwargs['label']
            del kwargs['label']

        f = open(os.path.join(self.directory, self.label + '.inpt'), 'w')

        # this finds it's way into the kwargs and isn't needed
        if 'directory' in kwargs:
            del kwargs['directory']

        # deal with the equivalent arguments, this gets complicated
        args = list(kwargs.copy())
        for arg in args:
            if arg in equivalencies:
                if arg == 'h':
                    continue
                else:
                    if equivalencies[arg] in kwargs or equivalencies[arg].lower() in kwargs:
                        raise InputError('Both {} and {} were input into the sparc '
                                         'calculator. These are equivalent arguments. '
                                         'Please only enter one of these.'.format(arg,
                                                                                  equivalencies[arg]))
                    kwargs[equivalencies[arg]] = kwargs.pop(arg)

        # make all kwargs upper case
        for arg in list(kwargs):
            if arg.upper() in default_parameters:
                kwargs[arg.upper()] = kwargs.pop(arg)

        ############## Begin writing the file ####################
        f.write('# Input File Generated By SPARC ASE Calculator #\n')

        # deal with the finite differnce grid
        if 'H' in kwargs:
            kwargs['h'] = kwargs.pop('H')
            #return None
        mesh_args = [kwargs.get(a) for a in ['MESH_SPARCING', 'h', 'FD_GRID']]
        inputs_check = [a is not None for a in mesh_args]
        if inputs_check.count(True) > 1:
            raise CalculatorSetupError('You can only specify one of the '
                                       'following: `h`, `FD_GRID`, `MESH_SPACING`')
        if inputs_check.count(True) == 0:
            raise CalculatorSetupError('You must specify one of the '
                                       'following: `h`, `FD_GRID`, `MESH_SPACING`')
            
        if 'h' in kwargs:
            kwargs['MESH_SPACING'] = kwargs['h']
        elif 'FD_GRID' in kwargs:

            fd_grid = self.interpret_grid_input(atoms, **kwargs)

            if fd_grid is not None:
                f.write('FD_GRID: {} {} {}\n'.format(*fd_grid))

        # Deal with the unit cell

        # this just saves a few levels of logic below
        if 'LATVEC' not in kwargs:
            kwargs['LATVEC'] = None
        if 'CELL' not in kwargs:
            kwargs['CELL'] = None

        # Scold/warn the user about using the CELL/LATVEC arguments
        if kwargs['LATVEC'] is not None and kwargs['CELL'] is None:
            raise InputError(
                'If LATVEC is input, you also must provide the CELL argument')

        if kwargs['LATVEC'] is None and kwargs['CELL'] is not None:
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
                    raise InputError(
                        'The value of LATVEC must have 9 elements (3x3)')
                lattice = np.array(kwargs['LATVEC'])
                lattice = np.reshape(lattice, (3, 3))
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
                    lattice = np.reshape(lattice, (3, 3))

            # Decipher CELL input
            if type(kwargs['CELL']) not in [str, list, tuple]:
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
                f.write('  {}'.format(format(cell_param / Bohr, ' .15f')))
            atoms.center()
            f.write('\n')
        else:
            f.write('CELL:')
            for length in atoms.get_cell_lengths_and_angles()[:3]:
                f.write('  {}'.format(format(length / Bohr, ' .15f')))
            f.write('\nLATVEC:')
            for cell_vec in atoms.cell:
                lat_vec = cell_vec / np.linalg.norm(cell_vec)
                f.write('\n')
                for element in lat_vec:
                    f.write('{}  '.format(format(float(element), ' .15f')))
            f.write('\n')

        # Deal with boundary conditions
        if 'BOUNDARY_CONDITION' in kwargs and 'BC' in kwargs:
            raise CalculatorSetupError('both BC and BOUNDARY_CONDITION'
                                       ' were set please only enter one'
                                       ' of these arguments')
        if 'BOUNDARY_CONDITION' in kwargs:
            warnings.warn('the BOUNDARY_CONDITION arguement is depricated'
                          ' please use the BC argument instead')
        else:
            if 'BC' in kwargs:
                iter_pbc = False
                # BC can either be an iterable with 3 boolean elements
                # (i.e. [True, True, True]) or a string to be put directly
                # in (i.e. 'P P P')
                if type(kwargs['BC']) == str:
                    if len(kwargs['BC'].split()) != 3:
                        raise CalculatorSetupError('string arguments to BC '
                                                   'must be 3 spaced boundary'
                                                   'conditions'
                                                   ' (i.e. \'P P P\')')
                    f.write('BC: ' + kwargs['BC'])
                elif '__iter__' in dir(kwargs['BC']):
                    iter_pbc = True
                    if len(kwargs['BC']) != 3:
                        raise CalculatorSetupError('the BC argument must have 3 '
                                                   'boolean elements (i.e. '
                                                   '[True, True, True]')
                    pbc = kwargs['BC']
            else:  # if no boundary condition arugments are present, use the
                # pbc in the atoms object
                iter_pbc = True
                pbc = atoms.get_pbc()
            if iter_pbc:
                pbc_str = ''
                for condition in pbc:
                    if condition:
                        pbc_str += ' P'
                    elif not condition:
                        pbc_str += ' D'
                    else:
                        raise CalculatorSetupError('only booleans may be '
                                                   'entered for boundary '
                                                   'conditions')
                f.write('BC:{}\n'.format(pbc_str))

        # convert input to usable format
        kpt_grid = self.interpret_kpoint_input(atoms, **kwargs)
        f.write('KPOINT_GRID: {} {} {}\n'.format(*kpt_grid))

        # check the spin situation
        if [float(a) for a in atoms.get_initial_magnetic_moments()] != [0] * len(atoms):
            spin_warn = str('inital magentic moments were set on the'
                            ' provided atoms, bu the spin polarization '
                            'flag has not been set. To run spin polarized'
                            ' pass in the flag `SPIN_TYP = 1`')
            if 'SPIN_TYP' in kwargs:
                if int(kwargs['SPIN_TYP']) != 1:
                    warnings.warn(spin_warn)
            else:
                warnings.warn(spin_warn)

        # xc should be put in separately
        if 'xc' in kwargs:
            kwargs['EXCHANGE_CORRELATION'] = kwargs['xc']
            if kwargs['xc'] == 'PBE':
                kwargs['xc'] = 'GGA_PBE'
            if kwargs['xc'] == 'LDA':
                kwargs['xc'] = 'LDA_PW'

            f.write('EXCHANGE_CORRELATION: ' +
                    default_parameters['EXCHANGE_CORRELATION'] + '\n')
            xc = default_parameters['EXCHANGE_CORRELATION']
        else:
            if 'EXCHANGE_CORRELATION' not in kwargs:
                f.write('EXCHANGE_CORRELATION: ' +
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
                xc = default_parameters['EXCHANGE_CORRELATION']
            else:
                if kwargs['EXCHANGE_CORRELATION'] == 'PBE':
                    kwargs['EXCHANGE_CORRELATION'] = 'GGA_PBE'
                if kwargs['EXCHANGE_CORRELATION'] == 'LDA':
                    kwargs['EXCHANGE_CORRELATION'] = 'LDA_PW'

                f.write('EXCHANGE_CORRELATION: ' +  # Note the Miss-spelling
                        kwargs['EXCHANGE_CORRELATION'] + '\n')
                xc = kwargs['EXCHANGE_CORRELATION']

        # you generally want this flag on if you're doing relaxation
        if 'RELAX_FLAG' in kwargs:
            if int(kwargs['RELAX_FLAG']) == 1:
                if 'PRINT_RELAXOUT' not in kwargs:
                    f.write('PRINT_RELAXOUT: 1\n')
                elif 'PRINT_RELAXOUT' == 0:
                    warnings.warn('When PRINT_RELAXOUT is set to 0 '
                                  'there is no way to recover the ')

        # turn on the flags for these print settings by default
        # it generally doesn't make sense to have these off
        for flag in ['PRINT_FORCES', 'PRINT_ATOMS']:
            if flag not in kwargs:
                f.write('{}: {}\n'.format(flag, '1'))

        # convert the Tols and such to eV and angstrom units
        hartree_inputs = ['SMEARING',
                          # 'TOL_SCF',
                          # 'TOL_POISSON',
                          # 'TOL_LANCZOS','TOL_PSEUDOCHARGE',
                          'SCF_ENERGY_ACC']
        bohr_inputs = ['MESH_SPACING']

        if 'TOL_RELAX' in kwargs:  # this is Ha/Bohr
            kwargs['TOL_RELAX'] /= Hartree / Bohr
        for setting in hartree_inputs:
            if setting in kwargs:
                kwargs[setting] /= Hartree
        for setting in bohr_inputs:
            if setting in kwargs:
                kwargs[setting] /= Bohr

        # all non-required inputs
        for arg, value in kwargs.items():
            if arg in special_inputs:
                continue
            if arg in equivalencies:
                continue
            f.write('{}: {}\n'.format(arg, value))

        f.close()

        # make the atomic inputs (.ion) file
        kwargs['pseudo_dir'] = self.get_pseudopotential_directory(
            xc=xc, **kwargs)

        outpath = os.path.join(self.directory, self.label)
        if kwargs.get('copy_psp') is None:
            copy_psp = True
        else:
            copy_psp = kwargs.get('copy_psp')
        write_ion(open(outpath + '.ion', 'w'),
                  atoms, pseudo_dir=kwargs['pseudo_dir'],
                  scaled=scaled, copy_psp=copy_psp)

    def interpret_grid_input(self, atoms, **kwargs):
        """
        helper function to translate whatever the user put in into
        a grid that can be directly written to an input file
        """
        if 'H' in kwargs:
            kwargs['h'] = kwargs.pop('H')
            #return None
        mesh_args = [kwargs.get(a) for a in ['MESH_SPARCING', 'h', 'FD_GRID']]
        inputs_check = [a is not None for a in mesh_args]
        if inputs_check.count(True) > 1:
            raise CalculatorSetupError('You can only specify one of the '
                                       'following: `h`, `FD_GRID`, `MESH_SPACING`')
        elif inputs_check.count(True)  == 0:
            raise CalculatorSetupError('You must specify one of the '
                                       'following: `h`, `FD_GRID`, `MESH_SPACING`')

        #if 'h' in kwargs:
        #    fd_grid = self.h2gpts(h=kwargs['h'], cell_cv=atoms.cell)
        # read the FD_GRID argument
        elif 'FD_GRID' in kwargs:
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
        else:
            raise InputError('one of the following must be specified:'
                             ' MESH_SPACING, h, or FD_GRID')
        return fd_grid

    def interpret_kpoint_input(self, atoms, **kwargs):
        """
        helper function to figure out what the kpoints input is and return
        it as an iterable
        """
        if kwargs.get('KPOINT_GRID'):
            if len(kwargs['KPOINT_GRID']) == 3:
                for kpoint in kwargs['KPOINT_GRID']:
                    if type(kpoint) != int:
                        raise InputError('when KPOINT_GRID is entered as an'
                                         ' iterable, the values must be in the'
                                         ' integer type (i.e. (4,4,4))')
                kpt_grid = kwargs['KPOINT_GRID']
            elif type(kwargs['KPOINT_GRID']) == str:
                if len(kwargs['KPOINT_GRID'].split()) != 3:
                    raise InputError('when KPOINT_GRID is entered as a string, it'
                                     ' must have 3 elements separated by spaces '
                                     '(i.e. \'4 4 4\')')
                kpt_grid = [int(a) for a in kwargs['KPOINT_GRID'].split()]
            elif len(kwargs['KPOINT_GRID']) != 3 or type(kwargs['KPOINT_GRID']) is not str:
                raise InputError('KPOINT_GRID must be either a length 3 object'
                                 ' (i.e. (4,4,4)) or a string (i.e. \'4 4 4 \')')
        else:
            kpt_grid = (1, 1, 1)
        return kpt_grid

        # check the spin situation
        if [float(a) for a in atoms.get_initial_magnetic_moments()] != [0] * len(atoms):
            spin_warn = str('inital magentic moments were set on the'
                            ' provided atoms, bu the spin polarization '
                            'flag has not been set. To run spin polarized'
                            ' pass in the flag `SPIN_TYP = 2`')
            if 'SPIN_TYP' in kwargs:
                if int(kwargs['SPIN_TYP']) != 2:
                    warnings.warn(spin_warn)
            else:
                warnings.warn(spin_warn)

    def get_pseudopotential_directory(self, **kwargs):
        if 'pseudo_dir' in kwargs.keys():
            return kwargs['pseudo_dir']
        elif 'pseudo_dir' not in kwargs.keys() and 'SPARC_PSP_PATH' in os.environ:
            kwargs['pseudo_dir'] = os.environ['SPARC_PSP_PATH']
            return os.environ['SPARC_PSP_PATH']
        elif 'pseudo_dir' not in kwargs.keys():
            # find where the defaults are
            current_file = inspect.getfile(inspect.currentframe())
            package_directory = os.path.dirname(os.path.abspath(current_file))
            psps_path = os.path.join(package_directory, 'pseudos')

            if 'LDA' in kwargs['xc']:
                psps_path = os.path.join(psps_path, 'LDA_pseudos')
            elif 'PBE' in kwargs['xc']:
                psps_path = os.path.join(psps_path, 'PBE_pseudos')
            kwargs['pseudo_dir'] = psps_path
            warnings.warn('No `pseudo_dir` argument was passed in and no '
                          '$SPARC_PSP_PATH environment variable was set '
                          'default pseudopotentials are being used.')
            return psps_path

    def get_nstates(self, atoms=None, **kwargs):
        if kwargs == {}:
            kwargs = self.parameters
        if atoms is None:
            atoms = self.atoms
        kwargs = kwargs.copy()

        if kwargs.get('nstates'):
            return kwargs.get('nstates')
        elif kwargs.get('NSTATES'):
            return kwargs.get('NSTATES')
        if 'xc' not in kwargs:
            if kwargs.get('EXCHANGE_CORRELATION'):
                kwargs['xc'] = kwargs.get('EXCHANGE_CORRELATION')

        psp_path = self.get_pseudopotential_directory(**kwargs)

        syms = list(set(atoms.get_chemical_symbols()))
        syms_valence_dict = {}
        for sym in syms:
            # read the pseudopotential file to find the number of valence
            with open(os.path.join(psp_path, sym + '.pot')) as f:
                lines = f.readlines()[:100]  # just grab the start
            for line in lines:
                if 'zion' in line:
                    valence = float(line.split()[1])
                    syms_valence_dict[sym] = valence
                    break  # stop once we find it
        tot_valence = sum([syms_valence_dict[a]
                           for a in atoms.get_chemical_symbols()])

        # per Qimen this is the formula in SPARC
        nstates = np.floor(tot_valence/2) * 1.2 + 5
        nstates = round(nstates)

        return nstates

    def setup_parallel_env(self):
        """
        just sets up the environment to have roughly optimal parallelization
        environment variables.
        """
        os.environ['MV2_ENABLE_AFFINITY'] = '1'
        os.environ['MV2_CPU_BINDING_POLICY'] = 'bunch'

    def generate_command(self):
        """
        simple function to seek out the location of `sparc` and `mpirun`
        executables in the PATH environment variable and build a command
        able to run sparc

        Parameters:
            None

        returns:
            command: (str)
                A command to run sparc

        """
        # initialize command as an empty string
        command = ''
        import shutil
        exe = shutil.which('sparc')
        if exe is None:
            raise EnvironmentError('No `sparc` executable command could be'
                                   ' found in the PATH environment variable')
        if 'PBS_NODEFILE' in os.environ:
            with open(os.environ['PBS_NODEFILE']) as f:
                np = len(f.readlines())
            mpi = which('mpirun')
            if mpi:
                command += '{} -np {}'.format(mpi, np)
        command += '{} -name PREFIX'.format(exe)

        return command

    def estimate_memory(self, atoms=None, units='GB', **kwargs):
        """
        a function to estimate the amount of memory required to run
        the selected calculation. This function takes in **kwargs,
        but if none are passed in, it will fall back on the parameters
        input when the class was instantiated
        """
        conversion_dict = {'MB': 1e-6, 'GB': 1e-9, 'B': 1, 'byte': 1,
                           'KB': 1e-3}
        if kwargs == {}:
            kwargs = self.parameters
        if atoms is None:
            atoms = self.atoms

        nstates = kwargs.get('NSTATES')
        if nstates is None:
            nstates = self.get_nstates(atoms=atoms, **kwargs)

        # some annoying code to figure out if it's a spin system
        spin_polarized = kwargs.get('nstates')
        if spin_polarized is not None:
            spin_polarized = int(spin_polarized)
        else:
            spin_polarized = 1
        if spin_polarized == 2:
            spin_factor = 2
        else:
            spin_factor = 1

        if 'MESH_SPACING' in kwargs:
            kwargs['h'] = kwargs.pop('MESH_SPACING')
        npoints = np.product(self.interpret_grid_input(atoms, **kwargs))

        kpt_grid = self.interpret_kpoint_input(atoms, **kwargs)
        kpt_factor = np.ceil(np.product(kpt_grid)/2)

        # this is a pretty generous over-estimate
        # TODO: check this function is working
        estimate = 5 * npoints * nstates * kpt_factor * spin_factor * 8  # bytes
        converted_estimate = estimate * conversion_dict[units]
        return converted_estimate

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes, command=None,
                  parallel=True, guess_command=False):
        # mash together all the output files for convienience
        self.concatinate_output()

        # sort out the command
        if command is None:
            if hasattr(self, 'command'):
                command = self.command
            elif 'ASE_SPARC_COMMAND' in os.environ:
                command = os.environ['ASE_SPARC_COMMAND']
                self.command = command
            else:
                pass
        else:
            self.command = command

        if 'forces' in properties:
            self.parameters['PRINT_FORCES'] = 1

        if 'stress' in properties:
            self.parameters['CALC_STRESS'] = 1

        if parallel == True:
            self.setup_parallel_env()
        if guess_command:
            self.command = self.generate_command()
        # more or less a copy of the parent class with small tweaks
        if atoms is not None:
            self.atoms = atoms

        self.write_input(self.atoms, **self.parameters)
        if self.command is None:
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)
        errorcode = subprocess.call(command, shell=True, cwd=self.directory)

        if errorcode:  # it will sometimes return an errorcode even though it finished
            if self.terminated_normally():
                pass
            else:
                path = os.path.abspath(self.directory)
                raise CalculationFailed('{} in {} returned an error: {}'
                                        .format(self.name, path, errorcode))
        self.read_results()

    def read_results(self):
        self.concatinate_output()
        # quick checks on the run to see if it finished
        if not self.terminated_normally():
            raise CalculationFailed('SPARC did not terminate normally'
                                    ' in the last run, please check the'
                                    ' output files to see what the error '
                                    'was')
        if not self.scf_converged():
            raise SCFError('The last SCF cycle of your run did not converge')
        parse_traj = False
        # regenerate the capitalized version of the input parameters
        kwargs = self.parameters
        # for arg in list(kwargs):
        #    if arg.upper() in default_parameters:
        #        kwargs[arg.upper()] = kwargs.pop(arg)

        # if it's MD or relaxation, we have to parse those files
        if 'MD_FLAG' in kwargs:
            if int(kwargs['MD_FLAG']) != 0 and kwargs['MD_FLAG'] is not None:
                parse_traj = True
        elif 'md_flag' in kwargs:
            if int('md_flag' in kwargs) != 0 and kwargs['md_flag'] is not None:
                parse_traj = True
        if 'RELAX_FLAG' in kwargs:
            if int(kwargs['RELAX_FLAG']) in [1, 2, 3]:
                parse_traj = True
        if 'relax_flag' in kwargs:
            if int(kwargs['relax_flag']) in [1, 2, 3]:
                parse_traj = True

        bundle = self.parse_output(parse_traj=parse_traj)
        if parse_traj:
            # update the atoms after relaxation/MD
            self.atoms.positions = bundle[0].get_positions()
            self.atoms.set_chemical_symbols(bundle[0].get_chemical_symbols())

    def terminated_normally(self):
        """
        reads the last few lines of a file to check if sparc finished
        normally
        """
        # check that the output file was made
        if not os.path.isfile(self.label + '.out'):
            return False
        # check that the termination block is there
        with open(self.label + '.out', 'r') as f:
            last_few_lines = f.readlines()[-10:]  # just to narrow the search
            for line in last_few_lines:
                if 'Acknowledgements: U.S. DOE (DE-SC0019410)' in line:
                    return True
        return False

    def scf_converged(self):
        """
        checks if the SCF converged in the last step of the .out file
        if the last step converged, but some failed in the middle of
        the run, it warns the user.
        """
        with open(self.label + '.out', 'r') as f:
            txt = f.read()
        txt = txt.split('SPARC')[-1]  # get the last run
        if 'did not converge to desired accuracy' in txt:
            some_failed = True
        else:
            return True
        last_step = txt.split('=' * 68)[-3]
        if 'did not converge to desired accuracy' in last_step:
            self.converged = False
            return False
        else:
            if some_failed:
                warnings.warn('At least one SCF cycle did not converge.'
                              ' Please check your results, as they may '
                              'be suspect')
                return True

    def get_scf_steps(self, include_uncompleted_last_step=False):
        """
        Gets the number of SCF steps in a geometric step

        inputs:

        include_uncompleted_last_step (bool):
            If set to True, the parser will count any uncompleted SCF cycles
            at the end of the file. This is only relevant if SPARC did not 
            terminate normally.

        returns:

        steps (list):
            A list containing the number of SCF steps for each geometric
            step in order
        """
        if self.results == {}:  # Check that SPARC has run
            return None
        f = open(self.label+'.out', 'r')
        out = f.read()
        # isolate the most recent run
        out = out.split('Input parameters')[-1]
        # cut the last run into individual SCF steps
        out_steps = out.split('SCF#')
        del out  # free up memory
        steps = []
        for step in out_steps:
            for line in step.split('\n'):
                # The total number of steps is printed at the end of each cycle
                if 'Total number of SCF' in line:
                    steps.append(int(line.split(':')[1]))
                    break
        return steps

    def get_geometric_steps(self, include_uncompleted_last_step=False):
        """
        Gets the number of geometric steps the run had

        Parameters:
            include_uncompleted_last_step (bool):
                If set to True, the parser will count any uncompleted SCF cycles
                at the end of the file. This is only relevant if SPARC did not 
                terminate normally.

        returns:
            num_steps (int):
                The number of geometric steps
        """
        inc = include_uncompleted_last_step  # because pep8
        steps = self.get_scf_steps(include_uncompleted_last_step=inc)
        num_steps = len(steps)
        return num_steps

    def get_runtime(self):
        """
        Parses the walltime from the SPARC output file
        """
        if self.results == {}:  # Check that SPARC has run
            return None
        with open(self.label + '.out', 'r') as f:
            txt = f.read()
        txt = txt.split('SPARC')[-1]  # get the most recent run
        walltime = re.findall('^.*Total walltime.*$', txt, re.MULTILINE)[-1]
        time = self.read_line(walltime, strip_text=True)
        return time

    def get_fermi_level(self):
        if hasattr(self, 'fermi_level'):
            return self.fermi_level
        else:
            self.calculate()

    def concatinate_output(self):
        """
        Combines together all outputs in the current directory.
        """
        files = os.listdir(self.directory)
        files.sort()
        for sufix in ['.out', '.relax', '.restart', '.aimd', '.geopt', '.static']:
            for item in files:
                if item.startswith(self.label) and sufix in item and \
                        item != self.label + sufix:
                    f = open(self.label + sufix, 'a')
                    g = open(os.path.join(self.directory, item), 'r')
                    text = g.read()
                    f.write('\n' + text)
                    os.remove(item)

    def read_line(self, text, typ=float, strip_text=False):
        """
        helper function to parse lines that are simple a colon separated
        keyword and value pair.

        Parameters:
            text (str):
                The line of text
            typ (type object):
                the type of the thing you'd like to parse
            strip_text (bool):
                set this to True if there is superfluous text at the
                end of the value

        returns:
            value:
                The value of the filed
        """
        if strip_text:
            if typ == float:
                return float(text.split()[-2])
            elif typ == int:
                return int(text.split()[-2])
            elif typ == str:
                return str(text.split(':')[-2])

        if typ == float:
            return float(text.split(':')[1])
        elif typ == int:
            return int(text.split(':')[1])
        elif typ == str:
            return str(text.split(':')[1])

    def parse_output(self, label=None, parse_traj=False,
                     return_results=False):
        """
        This function attempts to parse the output files of SPARC.
        You may ask it to parse the auxiliary atomic positions and
        forces files produced by relaxation or MD(.aimd and .geopt)
        This will use the files in self.label named files to reconstruct
        the inputs

        Parameters:
            parse_traj (bool):
                if set to True, it will attempt to parse the .aimd
                or .geopt files depening on if the .inpt file has
                the correct flags corresponding to MD or relaxation
                respectively

        returns:
           bundle (tuple)
                a tuple containing the information requested based on
                the values of `recover_input` and `parse_traj`

        """
        if label is not None:
            self.label = label
        with open(self.label + '.out', 'r') as f:
            text = f.read()
        text = text.split('SPARC')[-1]  # Get only the last run in the file
        text = text.split('*' * 75)  # break it into blocks

        """
        This is a guide to the different blocks in the `text` variable
        * the first two entires are the header
        * the third entry is input parameters
        * the fourth and fifth entries are the parallelization settings
        * the sixth is the header for the SCF runs
        * the seventh is the actual calculation
        """
        # get the details of the calculation
        self.results = {}
        step_split = text[-5].split('=' * 68)
        last_step = step_split[-1].splitlines()
        initialization = step_split[0].splitlines()
        atom_types = []
        n_atoms_of_type = []
        atom_valences = []
        atom_coordinates = []
        fractional = False
        forces = []
        stress = []

        for i, line in enumerate(initialization):
            # if 'oordinates of atoms' in line:
            #    if 'Fractional' in line:
            #        Fractional = True
            #    for j in range(i+1, i + n_atoms_of_type[-1] + 1):
            #        atom_coordinates.append(initialization[j].split())
            if 'Atom type' in line:
                atom_types.append(line.split()[-2])
                atom_valences.append(line.split()[-1])
            elif 'Pseudopotential' in line:
                pseudo = self.read_line(line, typ=str)
            elif 'Total number of electrons' in line:
                self.nvalence = self.read_line(line, typ=int)
            elif 'Total number of atoms' in line:
                natoms = self.read_line(line, typ=int)
            elif 'Number of atoms of type' in line:
                n_atoms_of_type.append(self.read_line(line, typ=int))

        # get the energy and forces from the last step
        for i, line in enumerate(last_step):
            if 'Pressure' in line:
                self.pressure = self.read_line(line, strip_text=True)
                self.pressure *= GPa
            elif 'Stress' in line:
                for j in range(i + 1, i + 4):
                    stress.append(last_step[j].split())
                stress = np.array(stress, dtype='float64')
                stress *= GPa
                assert np.shape(stress) == (3, 3)
                # convert to Voigt form
                stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                   stress[1, 2], stress[0, 2], stress[0, 1]])
                self.results['stress'] = stress
            elif 'Fermi level' in line:
                self.fermi_level = self.read_line(
                    line, strip_text=True) * Hartree
            elif 'Total free energy' in line:
                energy = self.read_line(line, strip_text=True) * Hartree
                self.results['energy'] = energy
                self.results['free_energy'] = energy
            elif 'Atomic forces' in line:
                for j in range(i + 1, i + natoms + 1):
                    forces.append(last_step[j].split())
                forces = np.array(forces, dtype='float64')
                forces *= Hartree / Bohr
                self.results['forces'] = forces
            elif 'Entropy*kb*T' in line:
                TS = self.read_line(line, strip_text=True) * Hartree
                # TODO: figure out what to do with this
                #self.results['energy'] -= TS * 0.5

        # Recover the original input parameters set,
        input_dict = self.parse_input_args(text[2])

        # parse the static file if it's an SCF calculation
        if 'MD_FLAG' not in input_dict and 'RELAX_FLAG' not in input_dict:
            with open(self.label + '.static') as f:
                static = f.read()
            static = static.split('*' * 75)[-1]
            atom_positions, static = static.split('\nTotal free energy (Ha): ')

            # TODO: think of a better way than this to split off the energy
            static = static.splitlines()
            total_energy = static[0]
            static = '\n'.join(static[1:])
            for line in atom_positions.splitlines()[1:]:
                if 'oordinates of' in line:
                    if 'Fractional' in line:
                        Fractional = True
                    continue
                atom_coordinates.append(line.split())
            if 'Atomic forces' in static:
                _, forces = static.split('Atomic forces (Ha/Bohr):\n')
                forces = forces.splitlines()[:len(atom_coordinates)]
                forces = [a.split() for a in forces]
                forces = np.array(forces, dtype='float64')
                forces *= Hartree / Bohr
                self.results['forces'] = forces

            if 'Stress ' in static:
                _, stress = static.split('Stress (GPa):')
                stress = stress.splitlines()[1:4]
                stress = np.array([a.split() for a in stress], dtype='float64')
                stress *= GPa
                assert np.shape(stress) == (3, 3)
                # convert to Voigt form
                stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                   stress[1, 2], stress[0, 2], stress[0, 1]])
                self.results['stress'] = stress
        if parse_traj:
            read_atoms = read_ion(open(self.label + '.ion', 'r'),
                                  recover_indices=False)
            inds = self.recover_index_order_from_ion_file(self.label)
            # you need a scrambled version of the indices to read the MD/relax files
            if 'MD_FLAG' in text[2]:
                # this function automatically unscrambles indices
                atoms = self.parse_MD(label=self.label, write_traj=True,
                                      pbc=read_atoms.pbc,
                                      cell=read_atoms.cell,
                                      chemical_symbols=read_atoms.get_chemical_symbols(),
                                      constraints=read_atoms.constraints,
                                      reorder=True)
                self.results['forces'] = atoms.get_forces()
                if 'stress' in atoms.calc.results:
                    self.results['stress'] = atoms.get_stress()
            elif 'RELAX_FLAG: 1' in text[2] or 'RELAX_FLAG: 3' in text[2]:
                # this function automatically unscrambles indices
                atoms = self.parse_relax(label=self.label, write_traj=True,
                                         pbc=read_atoms.pbc,
                                         cell=read_atoms.cell,
                                         chemical_symbols=read_atoms.get_chemical_symbols(),
                                         constraints=read_atoms.constraints,
                                         reorder=True)
                self.results['forces'] = atoms.get_forces()
                if 'stress' in atoms.calc.results:
                    self.results['stress'] = atoms.get_stress()
            else:
                # quickly unscramble
                atoms = read_atoms
                if len(inds) == len(atoms):
                    new_atoms = Atoms(['X'] * len(atoms),
                                      positions=[(0, 0, 0)] * len(atoms))
                    new_atoms.set_cell(atoms.cell)
                    if 'forces' in self.results:
                        forces = np.empty((len(atoms), 3))
                    # reassign indicies
                    for old_index, new_index in enumerate(inds):
                        new_atoms[new_index].symbol = atoms[old_index].symbol
                        new_atoms[new_index].position = atoms[old_index].position
                        new_atoms.pbc = atoms.pbc
                        if 'forces' in self.results:
                            forces[new_index] = self.results['forces'][old_index]
                    # make sure the atoms are correct
                    assert new_atoms.get_chemical_formula() == atoms.get_chemical_formula()
                    atoms = new_atoms
                    with open(self.label + '.static') as f:
                        static = f.read()
                    static = static.split('*' * 75)[-1]

                    if 'forces' in self.results:
                        self.results['forces'] = forces
                else:
                    warnings.warn('The SPARC ASE calculator was unable to '
                                  'reconstruct the atoms object from the input'
                                  ' files, this likely means that the .ion file'
                                  ' was modified during the run. Proceed with'
                                  ' Caution')

        bundle = []
        if parse_traj:
            bundle.append(atoms)
        bundle.append(input_dict)
        if return_results:
            bundle.append(self.results)
        return tuple(bundle)

    def parse_relax(self, label, write_traj=False,
                    pbc=False, cell=None,
                    chemical_symbols=[],
                    constraints=None,
                    reorder=True):
        if os.path.isfile(label + '.geopt'):
            f = open(label + '.geopt')
        elif os.path.isfile(label + '.restart'):
            f = open(label + '.restart')
        else:
            raise Exception('no .geopt or .restart file was found, make'
                            ' sure that you are running a relaxtion and'
                            ' that these files are not being deleted.')
        text = f.read()
        f.close()
        if cell is None and chemical_symbols == []:
            try:
                warnings.warn('attempting to rebuild atoms from'
                              'input file')
                atoms = read_ion(open(label + '.ion', 'r'))
                cell = atoms.cell
                chemical_symbols = atoms.get_chemical_symbols()
                pbc = atoms.pbc
            except:
                raise Exception('a SPARC .relax file cannot fully '
                                'define a system, you must either '
                                'input the chemical symbols and cell'
                                'or have a .ion and .inpt file for the'
                                ' system')
        # Parse out the steps
        if text == '':
            return None
        steps = text.split(':RELAXSTEP: 1')[-1]  # cut previous runs
        steps = steps.split(':RELAXSTEP:')
        if write_traj == False:
            steps = [steps[-1]]
        else:
            traj = Trajectory(label + '.traj', mode='w')

        if f.name.endswith('.restart'):
            steps = [steps[-1]]

        # Parse out the energies
        n_geometric = len(steps)
        with open(label + '.out', 'r') as f:
            s = f.read()
        s = s.split('SPARC')[-1]  # I think this makes searching faster
        s = re.findall('^.*Total free energy.*$', s, re.MULTILINE)
        engs = s[-n_geometric:]
        engs = [float(a.split()[-2]) * Hartree for a in engs]

        # build a traj file out of the steps
        for j, step in enumerate(steps):
            # Sometimes if it fails, a final empty step will exist
            if len(step.splitlines()) == 1:
                continue
            colons = step.split(':')
            positions = step.split(':')[4].strip().split('\n')
            forces = step.split(':')[6].strip().split('\n')
            frc = np.empty((len(forces), 3))
            stress = np.empty((3, 3))
            atoms = Atoms()
            for i, f in enumerate(forces):
                frc[i, :] = [float(a) * Hartree / Bohr for a in f.split()]
                atoms += Atom(chemical_symbols[i],
                              [float(a) * Bohr for a in positions[i].split()])
            stresses_found = False
            if 'STRESS' in step:
                stresses_found = True
                stress_index = colons.index('STRESS') + 1
                for i, s in enumerate(colons[stress_index].strip().split('\n')):
                    stress[i, :] = [float(a) * GPa for a in s.split()]
                # Convert of Voigt form
                stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                   stress[1, 2], stress[0, 2], stress[0, 1]])
            if reorder:
                inds = self.recover_index_order_from_ion_file(self.label)
                new_atoms = Atoms(['X'] * len(atoms))
                for old_index, new_index in enumerate(inds):
                    new_atoms[new_index].symbol = atoms[old_index].symbol
                    new_atoms[new_index].position = atoms[old_index].position
                # reorder the constraints
                if constraints is not None:
                    new_constraints = []
                    new_indices = []
                    for constraint in constraints:
                        if type(constraint).__name__ == 'FixAtoms':
                            for index in constraint.index:
                                new_indices.append(inds.index(index))
                            cons = constraint.copy()
                            new_constraints.append(cons)
                        elif type(constraint).__name__ in ['FixedLine', 'FixedPlan']:
                            for index in [constraint.a]:
                                new_index = inds[index]
                                cons = constraint.copy()
                                cons.index_shuffle(
                                    new_atoms, (index, new_index))
                                new_constraints.append(cons)
                        else:
                            warnings.warn('A constraint that is not of the types '
                                          'FixAtoms, FixedPlane, or FixedLine was '
                                          'found. These are not supported by SPARC'
                                          ' and will be ignored.')
                atoms = new_atoms
                atoms.set_constraint(new_constraints)
            elif constraints is not None:
                atoms.set_constraint(constraints)
            atoms.set_pbc(pbc)
            atoms.cell = cell
            if stresses_found:
                atoms.set_calculator(SinglePointCalculator(atoms, energy=engs[j],
                                                           forces=frc, stress=stress))
            else:
                atoms.set_calculator(SinglePointCalculator(atoms, energy=engs[j],
                                                           forces=frc))
            if write_traj == True:
                traj.write(atoms)
        return atoms

    def parse_MD(self, label, write_traj=False,
                 pbc=False, cell=None,
                 chemical_symbols=[],
                 constraints=None,
                 reorder=True):
        with open(label + '.aimd') as f:
            text = f.read()
        if cell is None and chemical_symbols == []:
            try:
                warnings.warn('attempting to rebuild atoms from'
                              'input file')
                atoms = read_ion(open(label + '.ion', 'r'))
                cell = atoms.cell
                chemical_symbols = atoms.get_chemical_symbols()
                pbc = atoms.pbc
            except:
                raise Exception('a SPARC .aimd file cannot fully '
                                'define a system, you must either '
                                'input the chemical symbols and cell'
                                'or have a .ion and .inpt file for the'
                                ' system')

        if text == '':
            return None
        # Parse out the steps
        steps = text.split(':MDSTEP: 1')[-1]  # get only the most recent run
        steps = steps.split(':MDSTEP:')
        if write_traj == False:
            steps = [steps[-1]]
        else:
            traj = Trajectory(label + '.traj', mode='w')

        # Parse out the energies
        n_images = len(steps)
        with open(label + '.out', 'r') as f:
            s = f.read()
        s = s.split('SPARC')[-1]  # I think this makes searching faster
        # regex to search for this line
        s = re.findall('^.*Total free energy.*$', s, re.MULTILINE)
        engs = s[-n_images:]
        engs = [float(a.split()[-2]) * Hartree for a in engs]

        # build a traj file out of the steps
        for j, step in enumerate(steps):
            # Find Indicies
            colons = step.split(':')
            #pos_index = colons.index('R(Bohr)') + 1
            #frc_index = colons.index('F(Ha/Bohr)') + 1
            #vel_index = colons.index('V(Bohr/atu)') + 1
            pos_index = colons.index('R') + 1
            frc_index = colons.index('F') + 1
            vel_index = colons.index('V') + 1
            # Parse the text
            positions = colons[pos_index].strip().splitlines()
            forces = colons[frc_index].strip().splitlines()
            velocities = colons[vel_index].strip().splitlines()
            # Initialize the arrays
            frc = np.empty((len(forces), 3))
            vel = np.empty((len(velocities), 3))
            stress = np.zeros((3, 3))
            atoms = Atoms()
            for i, f, v in zip(range(len(forces)), forces, velocities):
                frc[i, :] = [float(a) * Hartree / Bohr for a in f.split()]
                vel[i, :] = [float(a) * Bohr / fs for a in v.split()]
                atoms += Atom(chemical_symbols[i],
                              [float(a) * Bohr for a in positions[i].split()])
            stresses_found = False
            if 'STRESS' in step:
                stresses_found = True
                stress_index = colons.index('STRESS_TOT(GPa)') + 1
                for i, s in enumerate(colons[stress_index].strip().split('\n')):
                    stress[i, :] = [float(a) * GPa for a in s.split()]
                # Convert of Voigt form
                stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                   stress[1, 2], stress[0, 2], stress[0, 1]])

            if reorder:
                inds = self.recover_index_order_from_ion_file(self.label)
                new_atoms = Atoms(['X'] * len(atoms))
                for old_index, new_index in enumerate(inds):
                    new_atoms[new_index].symbol = atoms[old_index].symbol
                    new_atoms[new_index].position = atoms[old_index].position
                # reorder constraints
                if constraints is not None:
                    new_constraints = []
                    for constraint in constraints:
                        if type(constraint).__name__ == 'FixAtoms':
                            for index in constraint.index:
                                new_index = inds.index(index)
                                cons = constraint.copy()
                                cons.index_shuffle(
                                    new_atoms, [index, new_index])
                                new_constraints.append(cons)
                        elif type(constraint).__name__ in ['FixedLine', 'FixedPlan']:
                            for index in constraint.a:
                                new_index = inds[index]
                                cons = contraint.copy()
                                cons.index_shuffle(
                                    new_atoms, (index, new_index))
                                new_contraints.append(cons)
                        else:
                            warnings.warn('A constraint that is not of the types '
                                          'FixAtoms, FixedPlane, or FixedLine was '
                                          'found. These are not supported by SPARC'
                                          ' and will be ignored.')
                    constraints = new_constraints
                atoms = new_atoms

            atoms.set_velocities(vel)
            atoms.set_pbc(pbc)
            atoms.cell = cell
            if constraints is not None:
                atoms.set_constraint(constraints)
            if stresses_found:
                atoms.set_calculator(SinglePointCalculator(atoms,
                                                           energy=engs[j] *
                                                           len(atoms),
                                                           stress=stress,
                                                           forces=frc))
            else:
                atoms.set_calculator(SinglePointCalculator(atoms,
                                                           energy=engs[j] *
                                                           len(atoms),
                                                           forces=frc))

            if write_traj == True:
                traj.write(atoms)
        return atoms

    def parse_input_args(self, input_block):
        input_parameters = input_block.splitlines()
        input_dict = {}
        in_lattice_block = False
        for input_arg in input_parameters:
            # once we find LATVEC, we analyze the next 3 lines differently
            if 'LATVEC' in input_arg or in_lattice_block:
                if not 'block_line' in locals():
                    input_dict['LATVEC'] = []
                    in_lattice_block = True
                    block_line = 1
                    continue
                lat_vec = [float(a) for a in input_arg.strip().split()]
                input_dict['LATVEC'].append(lat_vec)
                if block_line == 3:
                    in_lattice_block = False
                    input_dict['LATVEC'] = np.array(input_dict['LATVEC'])
                    continue
                else:
                    block_line += 1
                    continue

            # This next conditional is just make sure this code doesn't
            # break if lines are added in the future that don't have colons
            if ':' not in input_arg or input_arg[0] == '#':
                continue

            kw, arg = input_arg.split(':')
            input_dict[kw.strip()] = arg.strip()
            if len(arg.split()) > 1:  # Some arguments are lists
                input_dict[kw.strip()] = arg.split()
        if input_dict.get('OUTPUT_FILE'):
            input_dict['label'] = input_dict['OUTPUT_FILE']
            del input_dict['OUTPUT_FILE']

        cell = [float(a) for a in input_dict['CELL']]
        cell = np.eye(3) * cell * Bohr
        # sort out the boundary conditions
        if 'BC' in input_dict:
            pbc = []
            for condition in input_dict['BC']:
                if condition == 'P':
                    pbc.append(True)
                elif condition == 'D':
                    pbc.append(False)
                elif condition not in ['P', 'D']:
                    raise Exception('a non-allowed boundary condition'
                                    ' was detected, please check your'
                                    ' input file and ensure only \'P\''
                                    ' or \'D\' were selected in the BC'
                                    ' arugument')
        # TODO: remove this when the arugument has been fully removed
        elif 'BOUNDARY_CONDITION' in input_dict:
            if input_dict['BOUNDARY_CONDITION'] == '2':
                pbc = [True, True, True]
            elif input_dict['BOUNDARY_CONDITION'] == '1':
                pbc = [False, False, False]
            elif input_dict['BOUNDARY_CONDITION'] == '3':
                pbc = [True, True, False]
            elif input_dict['BOUNDARY_CONDITION'] == '4':
                pbc = [True, False, False]
        return input_dict

    def recover_index_order_from_ion_file(self, label):
        """
        quickly parses the .ion file with the same label 
        to get the order of the indices so the MD or relaxtion
        file can be reordered
        """
        with open(label + '.ion', 'r') as f:
            txt = f.readlines()
        inds = []
        for line in txt:
            if '# index' in line:
                inds.append(int(line.split('# index')[1]))
        return inds

    @staticmethod
    def read(label):
        """
        Attempts to regenerate the SPARC calculator from a previous
        set of input/output files
        """
        calc = SPARC()
        #os.path.isfile(label + '.out')
        if os.path.isfile(label + '.out'):
            atoms, input_dict, results = calc.parse_output(label=label,
                                                           parse_traj=True,
                                                           return_results=True)
            if not calc.scf_converged():
                raise Exception('The SCF on this calculation did not converge')
        elif os.path.isfile(label + '.inpt') and os.path.isfile(label + '.ion'):
            atoms = read_ion(open(label + '.ion', 'r'))
            with open(label + '.inpt') as f:
                inputs = f.read()
            input_dict = calc.parse_input_args(inputs)
            results = {}
        else:
            raise Exception('no .out file was found')
        # just reinitialize now
        calc.__init__(atoms=atoms, **input_dict)
        calc.results = results
        return calc

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.reset()
        self.atoms = atoms

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def h2gpts(self, h, cell_cv, idiv=4):
        cell_cv = np.array(cell_cv)
        cell_lengths = np.linalg.norm(cell_cv, axis=1)
        grid = np.ceil(cell_lengths/h)
        grid = np.maximum(idiv, grid)
        return [int(a) for a in grid]

    def todict(self, only_nondefaults=False, return_atoms=True):
        """
        Coverts calculator object into a dictionary representation appropriate to be placed
        in a MongoDB. By default this returns only the settings that are not the default values.
        This function is based loosely off the todict function from the Kitchen Group's VASP
        environment. in modifying this function please attempt to conform to the naming 
        conventions found there:
        https://github.com/jkitchin/vasp/blob/master/vasp/vasp_core.py
        only_nondefaults: bool
        If set to True, only the non-default keyword arguments are returned. If False all
            key word arguements are returned
        """

        dict_version = OrderedDict()
        if self.results != {}:
            bundle = self.parse_output(parse_traj=True)
            atoms, input_dict = bundle
        else:
            input_dict = dict(**self.parameters)
        dict_version.update(**input_dict)

        # Check if SPARC terminated normally
        #dict_version['SPARC-completed'] = self.terminated_normally()
        #dict_version['path'] = os.path.abspath('.').split(os.sep)[1:]
        if self.results == {}:
            return dict_version
        for prop in self.implemented_properties:
            val = self.results.get(prop, None)
            dict_version[prop] = val
        f = self.results.get('forces', None)
        if f is not None:
            dict_version['fmax'] = max(np.abs(f.flatten()))
        s = self.results.get('stress', None)
        if s is not None:
            dict_version['smax'] = max(np.abs(s.flatten()))
        time = self.get_runtime()
        #if time is not None:
        #    dict_version['elapsed-time'] = time
        #steps = self.get_scf_steps()
        #if steps is not None:
        #    dict_version['SCF-steps'] = steps
        #if hasattr(self, 'nvalence'):
        #    dict_version['nvalence'] = self.nvalence
        #dict_version['name'] = 'SPARC-X'
        # Try to get a version. Since the code is in alpha, this is just the
        # time of the most recent commit.
        c_dir = os.getcwd()
        try:
            os.chdir(os.environ['ASE_SPARC_COMMAND'][:-5])  # rewrite later
            f = os.popen('git log | grep "Date"')
            recent_commit_date = f.readlines()[0]
            rc = recent_commit_date.split('Date:')[1]
            rc = rc.split('-')[0]
            dict_version['version'] = rc.strip()
            os.chdir(c_dir)
        except:
            os.chdir(c_dir)
        # rewrite
        for item in ['EXCHANGE_CORRELATION']:  # These should always be in
            try:
                dict_version[item] = self.parameters[item]
            except:
                pass
        for item in dict_version:
            if type(dict_version[item]) not in \
                [dict, list, str, float, int, None, tuple] and \
                    dict_version[item] is not None:  # converting numpy arrays to lists
                try:
                    dict_version[item] = dict_version[item].tolist()
                except:
                    pass

        if return_atoms == True:
            dict_version['atoms'] = self.atoms_dict(self.atoms)
        return dict_version

    def calc_to_mongo(self,
                      host='localhost',
                      port=27017,
                      database='atoms',
                      collection='atoms',
                      user=None,
                      password=None):
        """
        inserts a dictionary version of the calculator into a mongo database.

        """
        from .mongo import MongoDatabase, mongo_doc

        db = MongoDatabase(
            host=host,
            port=port,
            database=database,
            collection=collection,
            user=user,
            password=password)
        d = mongo_doc(self.get_atoms())
        id_ = db.write(d)
        print('entry was added to db with id {}'.format(id_))
        return id_

    def atoms_dict(self, atoms):
        d = OrderedDict(atoms=[{'symbol': atom.symbol,
                                'position': json.loads(encode(atom.position)),
                                'tag': atom.tag,
                                'index': atom.index,
                                'charge': atom.charge,
                                'momentum': json.loads(encode(atom.momentum)),
                                'magmom': atom.magmom}
                               for atom in atoms],
                        cell=atoms.cell,
                        pbc=atoms.pbc,
                        info=atoms.info,
                        constraints=[c.todict() for c in atoms.constraints])
        return d

    def dict_atoms(self, d):
        atoms = Atoms([Atom(atom['symbol'],
                            atom['position'],
                            tag=atom['tag'],
                            momentum=atom['momentum'],
                            magmom=atom['magmom'],
                            charge=atom['charge'])
                       for atom in d['atoms']],
                      cell=d['cell'],
                      pbc=d['pbc'],
                      info=d['info'],
                      constraint=[dict2constraint(c) for c in d['constraints']])
        return atoms
