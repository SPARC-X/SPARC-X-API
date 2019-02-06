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
required_manual = ['PSEUDOPOTENTIAL_FILE','CELL','EXCHAGE_CORRELATION','FD_GRID','NTYPES','PSEUDOPOTENTIAL_LOCAL','KPOINT_GRID']
default_parameters = {
            # 'label': 'sprc-calc',
            # 'calculation_directory':'sprk.log',

            'BOUNDARY_CONDITION': 2,
            'EXCHANGE_CORRELATION': 'LDA_PZ',  # 'LDA'
            'KPOINT_GRID': (1, 1, 1),
            'MIXING_PARAMETER': 0.30,
            'CHEN_DEGREE': 20,
            'NSTATES': None,
            'MAXIT_SCF': 100,
            'TOL_SCF': 1.00E-05,
            'TOL_POISSON': 1.00E-06,
            'TOL_LANCZOS': 1.00E-02,
            'TOL_PSEUDOCHARGE': 1.00E-08,
            'TWTIME': 999999999.000000,
            'MIXING_PARAMETER': 0.30,
            'MIXING_HISTORY': 7,
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
            'CHEB_DEGREE': 25,
            'NTYPES': None,

            'NP_KPOINT_PARAL': None,
            'NP_BAND_PARAL': None,
            'NP_DOMAIN_PARAL': None,
            'NP_DOMAIN_PHI_PARAL': None,

            'TOL_RELAX': 1.00E-03,
            'PRINT_RELAXOUT': 0,
            'RELAX_FLAG': 0,
            'RELAX_METHOD': None,
            'RELAX_MAXITER': 300,
            'NLCG_sigma': 0.500000,
            'L_HISTORY': 20,
            'L_FINIT_STP': 0.005000,
            'L_MAXMOV': 0.200000,
            'L_AUTOSCALE': 1,
            'L_LINEOPT': 1,
            'L_ICURV': 1.000000,

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
    """
    default_parameters = {
            # 'label': 'sprc-calc',
            # 'calculation_directory':'sprk.log',
            
            
            'BOUNDRY_CONDITION': 2,
            'EXCHANGE_CORRELATION': 'LDA_PZ',  # 'LDA'
            'KPOINT_GRID': (1, 1, 1),
            'MIXING_PARAMETER': 0.30,
            'CHEN_DEGREE': 20,
            'NSTATES': None,
            'MAXIT_SCF': 100,
            'TOL_SCF': 1.00E-05,
            'TOL_POISSON': 1.00E-08,
            'TOL_LANCZOS': 1.00E-02,
            'TOL_PSEUDOCHARGE': 1.00E-08,
            'MIXING_PARAMETER': 0.30,
            'MIXING_HISTORY': 7,
            'PULAY_FREQUENCY': 1,
            'PULAY_RESTART': 0,
            'REFERENCE_CUTOFF': 0.50,
            'RHO_TRIGER': 3,
            'PRINT_FORCES': 0,
            'PRINT_ATOMS': 0,
            'PRINT_EIGEN': 0,
            'PRINT_DENSITY': 0,
            'PRINT_RESTART_FQ': 1,
            'PSEUDOPOTENTIAL_LOCAL': 4,
            # 'PSEUDOPOTENTIAL_FILE': '../psdpots/psd_oncv_{}.pot',
            'OUTPUT_FILE': None,
            
            'NP_KPOINT_PARAL': None,
            'NP_BAND_PARAL': None,
            'NP_DOMAIN_PARAL': None,
            'NP_DOMAIN_PHI_PARAL': None,
            'CELL': None,
            'FD_GRID': None,
            'FD_ORDER': 12,

            'TOL_RELAX': 1.00E-03,
            'PRINT_RELAXOUT': 0,
            'RELAX_FLAG': 0,
            'RELAX_METHOD': 'NLCG',
            
                        }
    
    equivalencies = {
            'xc': 'EXCHANGE_CORRELATION',
            'kpts': 'KPOINT_GRID',
            'nbands': 'NSTATES',
            
            'gpts': 'FD_GRID'
            
            
            }
    misc = {'h': 'grid_spacing', }
   """ 
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
        self.num_calculations = 0
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
            
        self.set(**kwargs)

       # self.write_input(atoms=atoms,write_defaults=write_defaults, **kwargs)

    @staticmethod
    def write_input(atoms = None, write_defaults = False,
                    verbosity = 'normal', label = 'sprc-calc',
                    directory = '.',
                     **kwargs):
        ################ Check the input arguements ###################
        if atoms is None:
            raise Exception('Atoms object is required to create input file')
        #FileIOCalculator.write_input(self, atoms=atoms)

        if directory != os.curdir and not os.path.isdir(directory):
            os.makedirs(directory)

        os.chdir(directory)
        if [a == 90 for a in atoms.get_cell_lengths_and_angles()[3:]]\
               != [True, True, True]:  # check cell is a rectangle
           raise NotImplementedError("""Unit cells must be rectangular\
, non-orthogonoal cells are currently under development'""")

        # replace all arguments that are nNone by default with None in the input dict
        for key in default_parameters:  
            if key not in kwargs.keys() and default_parameters[key] == None:
                kwargs[key] = None

        f = open(os.path.join(directory, label) + '.inpt','w')
        
        ############## Begin writing the file ####################
        f.write('# Input File Generated By SPARC ASE Calculator #\n')
        
        ###### Requred Inputs Block
        # This section writes the minimal inputs needed to run SPARC

        f.write('NTYPES: '+ str(len(set(atoms.get_chemical_symbols()))) + '\n')
        # input the cell
        # if the system has no cell, give 6 A of space on each side
        #if 'CELL' not in kwargs:
        #    kwargs['CELL'] = None
        if kwargs['CELL'] == None:
            if round(float(np.linalg.norm(atoms.cell)), 1) == 0:
                f.write('CELL:')
                cell = np.eye(3) * (np.max(atoms.positions, axis=0) + (6, 6, 6))
                atoms.set_cell(cell)
                for cell_param in atoms.get_cell_lengths_and_angles()[0:3]:
                    f.write(' ' + str((cell_param) / Bohr))
                atoms.center()
            elif [a == 90 for a in atoms.get_cell_lengths_and_angles()[3:]] \
                   == [True, True, True]:  # check cell is a rectangle (again)
                f.write('CELL:')
                for length in atoms.get_cell_lengths_and_angles()[:3]:
                    f.write(' ' + str(length / Bohr))
            else:
                raise NotImplementedError("""Unit cells must be rectangular\
, non-orthogonoal cells are currently under development'""")
        else:
            if type(kwargs['CELL']) not in [str,list,tuple]:
                raise ValueError('CELL must be entered as a list, tuple, or space separted string')
            f.write('CELL:')
            if type(kwargs['CELL']) == str:
                #self.parameters['CELL'] = \
                #        [float(a) for a in kwargs['CELL'].split()]
                kwargs['CELL'] = kwargs['CELL'].split()
            kwargs['CELL'] = [float(a) for a in kwargs['CELL']]
            cell = kwargs['CELL']
            if len(cell) != 3:
                raise ValueError('The value of CELL must have 3 elements')
            for i,length in enumerate(cell):
                f.write(' ' + str(length))
            atoms.cell = np.eye(3) * np.array(cell) * Bohr
            """
        elif round(float(np.linalg.norm(atoms.cell)), 1) == 0:
            f.write('CELL:')
            cell = np.eye(3) * (np.max(atoms.positions, axis=0) + (6, 6, 6))
            atoms.set_cell(cell)
            for cell_param in atoms.get_cell_lengths_and_angles()[0:3]:
                f.write(' ' + str((cell_param) / Bohr))
            atoms.center()
        if [a == 90 for a in atoms.get_cell_lengths_and_angles()[3:]] \
               == [True, True, True]:  # check cell is a rectangle
            f.write('CELL:')
            for length in atoms.get_cell_lengths_and_angles()[:3]:
                f.write(' ' + str(length / Bohr))"""
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
        
        # deal with pseudopotential file path
        if kwargs['PSEUDOPOTENTIAL_FILE'] is not None:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for psp_path in kwargs['PSEUDOPOTENTIAL_FILE']:
                f.write(psp_path + ' ')
            f.write('\n')
            """
            If issue with pseudos not being in an absolute path and the
            limited length of file location is fixed this code will become
            useful
        elif 'PSP_PATH' in os.environ:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for element in sorted(list(set(atoms.get_chemical_symbols()))):
                psp_path = os.path.join(os.environ['PSP_PATH'],'psd_oncv_'+element+'.pot')
                f.write(psp_path+' ')
            f.write('\n')
            """
        elif 'PSP_PATH' in os.environ:
            f.write('PSEUDOPOTENTIAL_FILE: ')
            for element in sorted(list(set(atoms.get_chemical_symbols()))):
                pseudos_in_dir = [a for a in os.listdir(os.environ['PSP_PATH']) \
                            if a.endswith(element+'.pot')]
                filename = [a for a in os.listdir(os.environ['PSP_PATH']) \
                            if a.endswith(element+'.pot')][0]
                os.system('cp $PSP_PATH/' + filename + ' .')
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

        # xc should be put in separately
        if 'printDens' in os.environ['ASE_SPARC_COMMAND']:
            if 'EXCHANGE_CORRELATION' not in [a.upper() for a in kwargs]:
                f.write('EXCHAGE_CORRELATION: ' +  # Note the miss-spelling
                        default_parameters['EXCHANGE_CORRELATION'] + '\n')
            else:
                f.write('EXCHAGE_CORRELATION: ' +  # Note the Miss-spelling
                        kwargs['EXCHANGE_CORRELATION'] + '\n')
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
  
        ####### Non Minimal Input Arguements Block
        # This takes care of all the other parameters, this section does most of the work

        for key,value in kwargs.items():
            if key in equivalencies:  # handle ase style inputs such as 'xc'
                f.write(equivalencies[key]+ ': ' + str(value) + '\n')

            elif key in default_parameters.keys() and key not in required_manual:
                try: 
                    float(default_parameters[key])
                    if float(default_parameters[key]) == float(value):
                        continue
                except:
                    pass
                if (default_parameters[key.upper()] != value) \
                   and value != None:  # if it's not default, write it
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
        write_ion(open(os.path.join(directory, label) + '.ion','w'),atoms)

 
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
        if self.command is None:
            raise RuntimeError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        #This is designed only for PACE at the moment
        if 'PBS_NODEFILE' in os.environ:
            f = open(os.environ['PBS_NODEFILE'],'r')
            nodes = len(f.readlines())
            self.num_nodes = nodes
        #command = self.command.replace('PREFIX', self.prefix)
        #errorcode = subprocess.call(command, shell=True, cwd=self.directory)
            errorcode = 9999  # initialize a non-zero errorcode
            tries = 0
            while errorcode !=0:
            #time.sleep(2) # 2 second cushion on either side for safety, can be removed later
                errorcode = subprocess.call('mpirun '
                                    +'-env MV2_ENABLE_AFFINITY=1 -env '
                                    +'MV2_CPU_BINDING_POLICY=bunch'
                                    +' -np ' + str(self.num_nodes) + ' '
                                    +self.command + ' -name ' + self.prefix,
                                    #+' -log_summary > mpi.log'+
                                    shell=True, cwd=self.directory)
                self.concatinate_output()
                tries += 1
                if tries > 1000:
                    break
            #time.sleep(2)
        else:
            errorcode = subprocess.call(self.command + ' -name ' + self.prefix,
                                    shell=True, cwd=self.directory)

        if errorcode:
            raise RuntimeError('{} in {} returned an error: {}'
                               .format(self.name, self.directory, errorcode))
        self.concatinate_output()
        self.read_results()
        if self.converged == False:
            Warning('SPARC did not converge')
        self.num_calculations += 1
        
        
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
                self.atoms,params = parse_output(self.label)
                self.atoms.cell = cell
            
        output = self.label + '.out'
        # read and parse output
        f = open(output,'r')
        log_text = f.read()
        log_text = log_text.split('SPARC')[-1]  # isolate the most recent run
        body = log_text.rsplit('Timing info')[-2]

        # check that the SCF cycle converged
        conv_check = body.rsplit('Energy')[-2]
        conv_check = conv_check.split('\n')[-3]  # Parse the last step
        conv_check = float(conv_check.split()[-2])

        # get the convergence critera
        if 'TOL_SCF' in self.parameters.keys():
            conv_criteria = self.parameters['TOL_SCF']
        else:
            conv_criteria = default_parameters['TOL_SCF']
       
        # check convergence 
        if conv_check > float(conv_criteria):
            self.converged = False  
        else:
            self.converged = True
        
        # Parse out the total energy, energies are in Ha
        energy_force_block = body.rsplit('Energy')[-1]
        energy_force_block = energy_force_block.split('\n')
        output_energies_in_order = []
        #energy_force_block = body.rsplit('Energy and atomic forces')[-1]
        for energy in energy_force_block[2:8]:
            _,eng = energy.split(':')
            output_energies_in_order.append(float(eng.strip()[:-5]))
        
        # read forces, forces are in Ha/Bohr
        if 'PRINT_FORCES: 1' in log_text:
            forces = np.empty((len(self.atoms),3))
            for i,force in enumerate(energy_force_block[9:-5]):
                forces[i,:] = [float(a) for a in force.split()]

        # This is the order these values are in in the output        
        free_eng, band_struc_eng, xc_eng, self_corr_eng, \
        Entr_kbt, fermi_level = output_energies_in_order

        self.results = {
                'energy': free_eng * Hartree,
                'free_energy': free_eng * Hartree,
                }
        if 'PRINT_FORCES: 1' in log_text:
            self.results['forces'] =  forces * Hartree * Bohr
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
                if self.kwargs[item] == default_parameters[item] and only_nondefaults == True:
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
