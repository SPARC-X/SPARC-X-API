"""
author: Ben Comer (Georgia Tech)
"""
import os
import subprocess
from ase.calculators.calculator import FileIOCalculator
#from ase.utils.timing import Timer
from .utilities import h2gpts, parse_output
import numpy as np
from ase.units import Bohr, Hartree
from .ion import write_ion


all_properties = ['energy', 'forces', 'stress', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']


all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']
required_manual = ['PSEUDOPOTENTIAL_FILE','CELL','EXCHANGE_CORRELATION','FD_GRID','PSEUDOPOTENTIAL_LOCAL','KPOINT_GRID', 'LATVEC']
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

            'gpts': 'FD_GRID'


            }
misc = {'h': 'grid_spacing', 
    
        }


class SPARC(FileIOCalculator):
    """
    ase calculator class for SPARC. SPARC can be obtained from:
    https://github.com/SPARC-X/SPARC
    """
    
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, restart=None, ignore_bad_restart=False,
                 label='sprc-calc', atoms=None, command=None,
                 write_defaults=False, verbosity='normal', **kwargs):
        """
        restart (None):
            Not implemented
        ignore_bad_restart (None):
            Not implemented
        label (str):
            The prefix you would like sparc files to have in your directory
        atoms (ase atoms object):
            The atomic structure to perform a calculation on. If None, you 
            can still pass one in as an arguement 
            (i.e. SPARC().get_potential_energy([atoms object]).)
            alternatively you may attach the calculator to an atoms object
            (i.e. [atoms object].set_calculator(SPARC()))
        command (str):
            A command used to run SPARC-X. This can also be set with the environment
            variable ASE_SPARC_COMMAND.
        verbosity (str):
            The level of verbosity you would like from SPARC. options are
            'normal', 'low', and 'high'
        """
        if restart is not None:
            atoms, kwargs = parse_output(label)
            del kwargs['label']
        for key in kwargs:  # Make all keys upper case for convenience
            if key.upper() in default_parameters:
                kwargs[key.upper()] = kwargs.pop(key)
        for key in kwargs:
            if key not in list(default_parameters.keys()) + \
            list(equivalencies.keys()) + list(misc.keys()):# and \
            #key.upper() not in default_parameters.keys():
                raise TypeError('Unknown input parameter {}'.format(key))

        if 'RELAX_METHOD' in [a.upper() for a in kwargs.keys()]:
            if kwargs['RELAX_METHOD'] is not None:
                kwargs['RELAX_FLAG'] = 1
                kwargs['PRINT_RELAXOUT'] = 1
                if 'relax_flag' in kwargs.keys():
                    del kwargs['relax_flag']
                if 'print_relaxout' in kwargs.keys():
                    del kwargs['relax_flag']

        self.kwargs = kwargs  # Save the original user input
        #self.timer = Timer()
        self.initialized = False
        self.write_defaults = write_defaults
        self.verbosity = verbosity
        FileIOCalculator.__init__(self, restart, ignore_bad_restart, label,
                                  atoms, command)
        self.atoms = atoms
        if 'PRINT_FORCES' not in kwargs and 'PRINT_ATOMS' not in kwargs and \
        'PRINT_EIGEN' not in kwargs and 'PRINT_DENSITY' not in kwargs:
            if verbosity == 'low':
                pass
            else:
                kwargs['PRINT_FORCES'] = 1
                kwargs['PRINT_ATOMS'] = 1
                if verbosity == 'high':
                    kwargs['PRINT_EIGEN'] = 1
                    kwargs['PRINT_DENSITY'] = 1
        
        # convert smearing to hartree
        if 'SMEARING' in kwargs:
            if kwargs['SMEARING'] is not None:
                kwargs['SMEARING'] /= Hartree

        self.calc_to_database = self.calc_to_mongo
 
        self.set(**kwargs)

       # self.write_input(atoms=atoms,write_defaults=write_defaults, **kwargs)

    @staticmethod
    def write_input(atoms = None, write_defaults = False,
                    verbosity = 'normal', label = 'sprc-calc',
                    directory = '.', coordinates = 'fractional',
                     **kwargs):
        """
        A method to write a sparc input set. This includes the .inpt
        input file as well as the .ion positions file.

        parameters:
            atoms (ASE Atoms Object):
                An ase atoms object for the system of interest
            write_defaults (bool):
                If set to true, the code will write all the default
                values into the .inpt file
            verbosity (str):
                a flag to set the level of output for the SPARC code,
                options are 'low', 'normal', and 'high'. It is
                recommended this be left on 'normal'. Note, the low 
                verbosity does not print forces
            label (str):
                the prefix for the output filenames. The files will be
                named [label].inpt and [label].ion. Note that sparc's
                ouput file names will all begin with this prefix
            directory (str):
                The path of the directory where the output is to be
                saved
            coordinates (str):
                options are 'fractional' and 'absolute'. It is 
                recommended that users stick to fractional.
            
        """
        ################ Check the input arguements ###################
        if atoms is None:
            raise Exception('Atoms object is required to create input file')
        #FileIOCalculator.write_input(self, atoms=atoms)

        if directory != os.curdir and not os.path.isdir(directory):
            os.makedirs(directory)

        os.chdir(directory)

        # replace all arguments that are None by default with None in the input dict
        for key in default_parameters:  
            if key not in kwargs.keys() and default_parameters[key] == None:
                kwargs[key] = None

        f = open(os.path.join(directory, label) + '.inpt','w')
        
        ############## Begin writing the file ####################
        f.write('# Input File Generated By SPARC ASE Calculator #\n')
        
        ###### Requred Inputs Block
        # This section writes the minimal inputs needed to run SPARC

        # Scold/warn the user about using the CELL/LATVEC arguments
        if  kwargs['LATVEC'] is not None and  kwargs['CELL'] is None:
            raise ValueError('If LATVEC is input, you also must provide the CELL argument')

        if  kwargs['LATVEC'] is None and  kwargs['CELL'] is not None:
            UserWarning('The CELL argument was entered, but not the LATVEC argument. The unit cell in the input atoms object will be ignored and the cell will be assumed to be orthogonal')

        # Deal with CELL and LATVEC inputs
        if kwargs['CELL'] is not None:
            # If there's no LATVEC input, just assume it's an orthogonal unit cell
            if kwargs['LATVEC'] is None:
                    lattice = np.eye(3)
            else:
                # Decipher the LATVEC input
                try:  # try to deal with numpy arrays by making them lists
                    kwargs['LATVEC'] = list(kwargs['LATVEC'])
                except:
                    pass
                # Check that LATVEC is a native iterable
                if type(kwargs['LATVEC']) not in [str,list,tuple]:
                    raise ValueError('LATVEC must be entered as a list, tuple, or space and linebreak separated string')
                # Deal with space separated raw text entries
                elif type(kwargs['LATVEC']) == str:
                    kwargs['LATVEC'] = kwargs['LATVEC'].split()
                    if len(kwargs['LATVEC']) != 9:
                        raise ValueError('The value of LATVEC must have 9 elements (3x3)')
                    lattice = np.array(kwargs['LATVEC'])
                    lattice = np.reshape(lattice, (3,3))
                # Deal with the simple case of lists or tuples
                else:
                    lattice = np.array(kwargs['LATVEC'])
                    lattice = np.reshape(lattice,(3,3))
            
            # Decipher CELL input
            if type(kwargs['CELL']) not in [str,list,tuple]:
                raise ValueError('CELL must be entered as a list, tuple, or space separated string')
            if type(kwargs['CELL']) == str:
                kwargs['CELL'] = kwargs['CELL'].split()
            kwargs['CELL'] = [float(a) for a in kwargs['CELL']]
            cell = np.array(kwargs['CELL'])
            if len(cell) != 3:
                raise ValueError('The value of CELL must have 3 elements')
            atoms.cell = atoms.cell = (lattice.T * cell).T * Bohr


        # Use the cell of the atoms object to write the input           
        if round(float(np.linalg.norm(atoms.cell)), 1) == 0:
            f.write('CELL:')
            cell = np.eye(3) * (np.max(atoms.positions, axis=0) + (6, 6, 6))
            atoms.set_cell(cell)
            for cell_param in atoms.get_cell_lengths_and_angles()[0:3]:
                f.write(' ' + str((cell_param) / Bohr))
            atoms.center()
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
        # input finite difference grid
        f.write('FD_GRID:')

        if kwargs['FD_GRID'] is not None:  # fix 
            if type(kwargs['FD_GRID']) == str:
                if len(n_pts.split(' ')) != 3:
                    raise Exception('if FD_GRID is entered as a string, it must be 3 numbers separated by spaces (i.e. 30 30 30).')
                f.write(n_pts)
            else:
                for n_pts in kwargs['FD_GRID']:
                    f.write(' ' + str(n_pts))
        elif 'h' in kwargs:
            fd_grid = h2gpts(kwargs['h'], atoms.cell, idiv=1)
            Warning('when utilizing h, grid points are set to approximately match the value of h. Enetered h values will not yield grids with spacing of exactly h')
            for n_pts in fd_grid:
              f.write(' ' + str(n_pts))
        elif 'gpts' in kwargs:
            for n_pts in kwargs['gpts']:
              f.write(' ' + str(n_pts))
        else:  # default of 0.2 A was chosen
            if 'cell' not in locals():
                cell = atoms.cell
            fd_grid = h2gpts(0.2, cell, idiv=1)
            for n_pts in fd_grid:
               f.write(' ' + str(n_pts))
        f.write('\n')

        # kpoints
        if 'kpts' in kwargs:
            if 'KPOINT_GRID' in kwargs:
                Warning('both kpts and KPOINT_GRID were input, defaulting to kpts')
            kwargs['KPOINT_GRID'] = kwargs['kpts']
        if 'KPOINT_GRID' in kwargs:
            if kwargs['KPOINT_GRID'] is not None:
                f.write('KPOINT_GRID: ')
                if len(kwargs['KPOINT_GRID']) == 3:
                    for kpoint in kwargs['KPOINT_GRID']:
                        if type(kpoint) != int:
                            raise Exception('when KPOINT_GRID is entered as an iterable, the values must be in the integer type (i.e. (4,4,4))')
                        f.write(str(kpoint) + ' ')
                elif type(kwargs['KPOINT_GRID']) == str:
                    if len(kwargs['KPOINT_GRID'].split(' ')) != 3:
                        raise Exception('when KPOINT_GRID is entered as a string, it must have 3 elements separated by spaces (i.e. \'4 4 4\')')
                    f.write(kwargs['KPOINT_GRID'])
                elif len(kwargs['KPOINT_GRID']) != 3 or type(kwargs['KPOINT_GRID']) is not str:
                    raise Exception('KPOINT_GRID must be either a length 3 object (i.e. (4,4,4)) or a string (i.e. \'4 4 4 \')')
                f.write('\n')
        
        # deal with pseudopotential file path
        """
        if kwargs['PSEUDOPOTENTIAL_FILE'] is not None:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for psp_path in kwargs['PSEUDOPOTENTIAL_FILE']:
                f.write(psp_path + ' ')
            f.write('\n')
        """
        """
            If issue with pseudos not being in an absolute path and the
            limited length of file location is fixed this code will become
            useful
        elif 'SPARC_PSP_PATH' in os.environ:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for element in sorted(list(set(atoms.get_chemical_symbols()))):
                psp_path = os.path.join(os.environ['SPARC_PSP_PATH'],'psd_oncv_'+element+'.pot')
                f.write(psp_path+' ')
            f.write('\n')
        """
        """
        elif 'SPARC_PSP_PATH' in os.environ:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for element in sorted(list(set(atoms.get_chemical_symbols()))):
                pseudos_in_dir = [a for a in os.listdir(os.environ['SPARC_PSP_PATH']) \
                            if a.endswith(element+'.pot')]
                filename = [a for a in os.listdir(os.environ['SPARC_PSP_PATH']) \
                            if a.endswith(element+'.pot')][0]
                os.system('cp $SPARC_PSP_PATH/' + filename + ' .')
                if pseudos_in_dir != [] and len(pseudos_in_dir) > 1:
                    f.write(pseudos_in_dir[0] + ' ')
                else:
                    psp_path = filename
                    f.write(psp_path + ' ')
            f.write('\n')
        elif write_defaults is True:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for element in sorted(list(set(atoms.get_chemical_symbols()))):
                f.write('../psdpots/psd_oncv_' + element + '.pot ')
            f.write('\n')
        """
        # xc should be put in separately
        if 'xc' in kwargs:
            kwargs['EXCHANGE_CORRELATION'] = kwargs['xc']
            """
        if 'printDens' in os.environ['ASE_SPARC_COMMAND']:
            if 'EXCHANGE_CORRELATION' not in [a.upper() for a in kwargs]:
                f.write('EXCHAGE_CORRELATION: ' +  # Note the miss-spelling
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
            else:
                f.write('EXCHAGE_CORRELATION: ' +  # Note the Miss-spelling
                        kwargs['EXCHANGE_CORRELATION'] + '\n')
            """
        else:
            if 'EXCHANGE_CORRELATION' not in [a.upper() for a in kwargs]:
                f.write('EXCHANGE_CORRELATION: ' +  # Note the miss-spelling
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
            else:
                f.write('EXCHANGE_CORRELATION: ' +  # Note the Miss-spelling
                        kwargs['EXCHANGE_CORRELATION'] + '\n')

        # check if the user has defined some print settings, these override verbosity inputs
        non_standard_verbosity_input = [a for a in kwargs if 'PRINT' in a]
        non_standard_verbosity_input = non_standard_verbosity_input != []
        if verbosity == 'low' and not non_standard_verbosity_input:
            pass
        elif verbosity == 'normal' and not non_standard_verbosity_input:
            f.write('PRINT_FORCES: 1\nPRINT_ATOMS: 1\n')
        elif verbosity == 'high' and not non_standard_verbosity_input:
            f.write('PRINT_FORCES: 1\nPRINT_ATOMS: 1\n \
                    PRINT_EIGEN: 1\nPRINT_DENSITY: 1\n')
  
        ####### Non Minimal Input Arguments Block
        # This takes care of all the other parameters, this section does most of the work

        for key,value in kwargs.items():
            if key in equivalencies:  # handle ase style inputs such as 'xc'
                if equivalencies[key] in required_manual:
                    continue
                f.write(equivalencies[key]+ ': ' + str(value) + '\n')

            elif key in default_parameters.keys() and key not in required_manual:
                try: 
                    float(default_parameters[key])
                    if float(default_parameters[key]) == float(value):
                        continue
                except:
                    pass
                if (default_parameters[key.upper()] != value) \
                   and value is not None:  # if it's not default, write it
                    if type(value) in [tuple, list]:
                        f.write(key)
                        for i in value:
                            f.write(' ' + str(i))
                        f.write('\n')
                    else:
                        f.write(key.upper() + ': ' + str(value) + '\n')
            #f.write('\n')
        #If the user wants the defaults explicitly in the input file
        if write_defaults is True:
            for key, value in default_parameters.items():
                if key not in kwargs.keys() and key !='EXCHANGE_CORRELATION' \
                and value != None:
                    f.write(key + ': '+str(value))
                f.write('\n')
        #write the atomic positions file (.ion file)
        write_ion(open(os.path.join(directory, label) + '.ion','w'),
                  atoms, pseudo_dir = os.environ['SPARC_PSP_PATH'],
                  scaled = False)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        """Do the calculation.

        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        """
        #FileIOCalculator.calculate(self,atoms=atoms)
        if atoms is not None:
            self.atoms = atoms.copy()
        else:
            atoms=self.atoms
        self.write_input(atoms = atoms, verbosity = self.verbosity,
                         label = self.label,directory = self.directory, 
                         **self.parameters)
        if 'PBS_NP' in os.environ:
            if float(os.environ['PBS_NP']) > 1:
                os.environ['MV2_ENABLE_AFFINITY'] = '1'
                os.environ['MV2_CPU_BINDING_POLICY'] = 'bunch'
        if self.command is None:
            raise RuntimeError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        #errorcode = subprocess.call(command, shell=True, cwd=self.directory)
            errorcode = 20  # initialize a non-zero errorcode
            tries = 0
            # random codes failures at the end have been mitigated with
            # a while loop
            command = self.command.replace('PREFIX', self.prefix)

            while errorcode !=0:
                #errorcode = subprocess.call('mpirun '
                                    #+'-env MV2_ENABLE_AFFINITY=1 -env '
                                    #+'MV2_CPU_BINDING_POLICY=bunch'
                                    #+' -np ' + str(self.num_nodes) + ' '
                                    #+self.command + ' -name ' + self.prefix,
                                    #+' -log_summary > mpi.log'+
                                    #shell=True, cwd=self.directory)
                errorcode = subprocess.call(command, shell = True,
                                            cwd = self.directory)
                self.concatinate_output()
                tries += 1
                if tries > 4:
                    break
        else:
            command = self.command.replace('PREFIX', self.prefix)
            errorcode = 20  # initialize a non-zero errorcode
            tries = 0

            while errorcode != 0:
                errorcode = subprocess.call(command,
                                    shell=True, cwd=self.directory)
                self.concatinate_output()
                tries += 1
                if tries > 4:
                    break

        if errorcode:
            raise RuntimeError('{} in {} returned an error: {}'
                               .format(self.name, self.directory, errorcode))
        self.concatinate_output()
        self.read_results()
        if self.converged == False:
            Warning('SPARC did not converge')
        
    def read_results(self):
        """"
        # find the most recently modified output file to read
        time = -10**6
        for output_file in os.listdir('./'):
            if self.label in output_file and 'out' in output_file:
                new_time = os.path.getmtime(output_file)
                if new_time > time:
                    time = new_time
                    output = output_file
        """
        # If the relaxation built into SPARC was used, read the results
        if 'RELAX_FLAG' in self.parameters:
            if self.parameters['RELAX_FLAG'] == 1:
                cell = self.atoms.cell
                atoms, params = parse_output(self.label,
                                                  calc_type = 'relax')
                if atoms is not None:
                    self.atoms.cell = cell
        if 'MD_FLAG' in self.parameters:
            if self.parameters['MD_FLAG'] == 1:
                cell = self.atoms.cell
                atoms, params = parse_output(self.label,
                                                  calc_type = 'MD')
                if atoms is not None:
                    self.atoms.cell = cell

            
        output = self.label + '.out'

        # read and parse output
        f = open(output,'r')
        log_text = f.read()
        log_text = log_text.split('SPARC')[-1]  # isolate the most recent run
        body = log_text.rsplit('Timing info')[-2]

        # check that the SCF cycle converged
        conv_check = body.rsplit('Energy')[-2]
        if 'did not converge to desired accuracy' in conv_check:
            self.converged = False
        else:
            self.converged = True

        # Parse out the total energy, energies are in Ha
        energy_force_block = body.rsplit('Energy')[-1]
        energy_force_block = energy_force_block.split('\n')
        output_energies_in_order = []
        #energy_force_block = body.rsplit('Energy and atomic forces')[-1]
        for energy in energy_force_block[2:9]:
            _,eng = energy.split(':')
            eng,_ = eng.split('(')
            output_energies_in_order.append(float(eng.strip()))
        
        # read forces, forces are in Ha/Bohr
        if 'PRINT_FORCES: 1' in log_text:
            forces = np.empty((len(self.atoms),3))
            for i,force in enumerate(energy_force_block[10:-5]):
                forces[i,:] = [float(a) for a in force.split()]

        # This is the order these values are in in the output        
        p_atom_eng, free_eng, band_struc_eng, xc_eng, \
        self_corr_eng, Entr_kbt, fermi_level \
        = output_energies_in_order

        self.results = {
                'energy': free_eng * Hartree,
                'free_energy': free_eng * Hartree,
                }
        if 'PRINT_FORCES: 1' in log_text:
            self.results['forces'] =  forces * Hartree / Bohr
        #print('output read')
        self.fermi_level = fermi_level
       

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.reset()
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def get_fermi_level(self):
        if self.results == {}:
            self.calculate()
        return self.fermi_level

    def get_forces(self, atom=None):
        #print('force call')
        if self.parameters['PRINT_FORCES']==0:
            raise Exception('verbosity must be set to normal, or PRINT_FORCES must be set to 1 to call get_forces')
        if self.results == {}:
            self.calculate()
        return self.results['forces']
        #return self.get_property('forces', atoms)

    def concatinate_output(self):
        """
        Combines together all outputs in the current directory.
        """
        files = os.listdir('.')
        files.sort()
        for sufix in ['.out','.relax','.restart']:
            for item in files:
                if item.startswith(self.label) and sufix in item and \
                item != self.label + sufix:
                    f = open(self.label + sufix, 'a')
                    g = open(item, 'r')
                    text = g.read()
                    f.write('\n' + text)
                    os.remove(item)

    def get_runtime(self):
        """
        Parses the walltime from the SPARC output file
        """
        if self.results == {}:  # Check that SPARC has run
            return None
        f = os.popen('grep "Total walltime" ' + self.label + '.out')
        time = f.readlines()[-1].split(':')[1]
        time = float(time.split('sec')[0].strip())
        f.close()
        return time

    def get_geometric_steps(self):
        """
        Gets the number of geometric steps run in SPARC's internal optimizer

        returns:
            steps (int):
                The number of geometric steps run in SPARC's internal optimizer
        """
        if self.results == {}: # Check that SPARC has run
            return None
        if self.parameters['RELAX_FLAG'] == 1:
            f = os.popen('grep "SCF#" ' + self.label + '.out')
            steps = f.readlines()[-1].split('SCF#')[1]
            steps = int(steps.split(')'))
            return steps + 1
        else:
            return None            

    def get_scf_steps(self, step_num = -1, 
                      include_uncompleted_last_step = False):
        """
        Gets the number of SCF steps in a geometric step

        inputs:

        step_num (int):
            The step number for which the number of SCF steps will be returned
        include_uncompleted_last_step (bool):
            If set to True, the parser will count any uncompleted SCF cycles
            at the end of the file. This is only relevant if SPARC did not 
            terminate normally.

        returns:
        
        steps (int):
            The number of SCF steps in the chosen geometric step
        """
        if self.results == {}: # Check that SPARC has run
            return None
        f = open(self.label+'.out','r')
        out = f.read()
        out = out.split('SCF#')
        if step_num < 0:
            steps_block = out[step_num].split('=' * 68)[1]
        else:
            steps_block = out[step_num+1].split('=' * 68)[1]
        num_steps = len(steps_block.split('\n')) - 3
        return num_steps

    def todict(self, only_nondefaults = False, return_atoms = True):
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
    
        from collections import OrderedDict
        import getpass 
        dict_version = OrderedDict()
        for item in default_parameters:  # rewrite this section
            if item in self.kwargs.keys():
                
                if item == 'LATVEC':
                    dict_version[item] = self.kwargs[item]
                elif self.kwargs[item] == default_parameters[item] and\
                   only_nondefaults == True:
                    pass
                elif type(self.kwargs[item]) == dict:
                    dict_version[item] = OrderedDict(self.kwargs[item])
                else:
                    dict_version[item] = self.kwargs[item]
                
            elif only_nondefaults == False:
                dict_version[item] = default_parameters[item]
        if hasattr(self,'converged'):
            dict_version['SCF-converged'] = self.converged
        f = os.popen('tail ' + self.label + '.out')
        

        # Check if SPARC terminated normally
        if ' U.S. National Science' in f.readlines()[-2]:
            dict_version['SPARC-completed'] = True
        dict_version['path'] = os.path.abspath('.').split(os.sep)[1:]
        if self.results == {}:
            return dict_version
        for prop in self.implemented_properties:
            val = self.results.get(prop, None)
            dict_version[prop] = val 
        f = self.results.get('forces', None)
        if f is not None:
            dict_version['fmax'] = max(np.abs(f.flatten()))
        #s = self.results.get('stress', None)
        #if s is not None:
        #    dict_version['smax'] = max(np.abs(s.flatten()))
        time = self.get_runtime()
        if time is not None:
            dict_version['elapsed-time'] = time
        steps = self.get_scf_steps()
        if steps is not None:
            dict_version['SCF-steps'] = steps
        f = os.popen('grep "Total number of electrons" {}.out | tail -1'.format(self.label))
        if f is not None:
            txt = f.read()
            dict_version['nvalence'] = float(txt.split(':')[1].strip())
        dict_version['name'] = 'SPARC-X'
        # Try to get a version. Since the code is in alpha, this is just the
        # time of the most recent commit.
        c_dir = os.getcwd()
        try:
            os.chdir(os.environ['ASE_SPARC_COMMAND'][:-5]) #  rewrite later
            f = os.popen('git log | grep "Date"')
            recent_commit_date = f.readlines()[0]
            rc = recent_commit_date.split('Date:')[1]
            rc = rc.split('-')[0]
            dict_version['version'] = rc.strip()
            os.chdir(c_dir)
        except:
            os.chdir(c_dir)
        for item in ['EXCHANGE_CORRELATION','EXCHAGE_CORRELATION']: # These should always be in
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
            from .utilities import atoms_dict
            dict_version['atoms'] = atoms_dict(self.atoms)
        return dict_version

    @staticmethod
    def fromdict(atoms = None, **kwargs):
        """
        Regenerates the sparc calculator from a dictionary version
        """
    #   other_keys = [None, ignore_bad_restart=False,
    #                 label='sprc-calc', atoms=None, command=None,
    #                 write_defaults=False, verbosity='normal',]
        extraneous_args = [a for a in kwargs.keys() if a not in default_parameters.keys()]
        from .utilities import dict_atoms
        atoms = dict_atoms(atoms)
        for key in extraneous_args:
            del kwargs[key]
        return SPARC(atoms = atoms, **kwargs)
            

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

#   def get_number_of_grid_points(): implement in the future
