from .sparc import SPARC, default_parameters
from .utilities import atoms_dict, dict_atoms
import os
from fireworks import explicit_serialize, FiretaskBase, FWAction, Firework, Workflow
from ase.atoms import Atoms
from ase.atom import Atom
from .mongo import mongo_atoms_doc, mongo_doc_atoms, MongoDatabase
from pymatgen.io.ase import AseAtomsAdaptor
from collections import OrderedDict
from ase.io.jsonio import encode
import json


@explicit_serialize
class RunSparcASE(FiretaskBase):
    """
    Runs SPARC using ASE based on some input dictionary

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        
    """
    required_params = ['atoms', 'parameter_dict','to_db','identifier','db_file']
    optional_params = ['sparc_command']

    _fw_name = 'Run SPARC-X'
    def run_task(self,fw_spec):
        #if self.get("expand_vars"):
        #    os.environ['ASE_SPARC_COMMAND'] = fw_spec['sparc_command']
        #if fw_spec['psuedo_potentials_path'] is not None:
        #    os.environ['PSP_PATH'] = fw_spec['psuedo_potentials_path']
        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        calc = SPARC(**parameter_dict)
        atoms.set_calculator(calc)
        calc.calculate()
        if self['to_db'] == True:
            from json import load
            db_info = load(open(self['db_file'],'r'))
            id_ = calc.calc_to_mongo(**db_info )
            if self['identifier'] is not None:
                db = MongoDatabase(**db_info)
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)
        formula = atoms.get_chemical_formula()

        formula = atoms.get_chemical_formula()
        try:
            if self['db_file'] == None:
                os.system('mkdir /nv/hp13/bcomer3/data/density_outputs/' + self['identifier'])
                os.system('mv * /nv/hp13/bcomer3/data/density_outputs/' + self['identifier'] + '/')
        except:
            pass


class Sparc_SCF_FW(Firework):
    def __init__(self, atoms, 
                 parameters = default_parameters,  
                 name = 'SPARC SCF',
                 sparc_command = None,
                 psuedo_potentials_path = None,
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 **kwargs
                 ):  
        """
        Runs an SCF calculation in SPARC on the given structure

        Args:
            atoms (Atoms Object or Pymatgen Structure): Input structure
        
            name (str): Name for the Firework
        
            parameters (Dict): the input parameters for SPARC
        
            sparc_command (str): optional, a command used by ase to run SPARC.
                This can also be set with the $ASE_SPARC_COMMAND environment
                variable.

            psuedo_potentials_path (str): optional, a path to the pseudopotentials
                you'd like used. This can also be set using the $PSP_PATH environment
                Variable.

            db_file (str): Not implemented

            parents ([Firework]): Parents of this Firework
        """
    
        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(RunSparcASE(atoms = atoms, parameter_dict = parameters,
             sparc_command = sparc_command, 
             psuedo_potentials_path = psuedo_potentials_path,
             to_db = to_db,
             db_file = db_file,
             identifier = identifier))
        super(Sparc_SCF_FW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

def get_sparc_runs(structures, parameters = default_parameters,
                                sparc_command = None,
                                to_db = True,
                                db_file = None,
                                optimize = False,
                                psuedo_potentials_path = None,
                                identifiers = None):
    fws = []
    if type(parameters) != list:  # If no list of parameters is given, use the same for all
        parameters = [parameters] * len(structures) 
    if type(identifiers) != list:
        identifiers = [identifiers]
    if optimize == False:
        FW = Sparc_SCF_FW
    else:
        FW = Sparc_Optimize_FW 
    for struct, param, identifier in zip(structures, parameters, identifiers):
        name = struct.get_chemical_formula()
        fws.append(FW(atoms_dict(struct),param,
                    sparc_command = sparc_command,
                    psuedo_potentials_path = psuedo_potentials_path,
                    to_db = to_db,
                    db_file = db_file,
                    identifier = identifier))
    return Workflow(fws, name="{} tests wf, e.g.,".format(len(fws)))
    #return Workflow(fws, name="{} tests wf, e.g., {}".format(len(fws), fws[0].name))


@explicit_serialize
class OptimizeLattice(FiretaskBase):
    """
    Runs SPARC using ASE based on some input dictionary

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        
    """
    required_params = ['atoms', 'parameter_dict','to_db','identifier']
    optional_params = ['sparc_command']

    _fw_name = 'Write Sparc Input From Dict'
    def run_task(self,fw_spec):
        #if self.get("expand_vars"):
        #    os.environ['ASE_SPARC_COMMAND'] = fw_spec['sparc_command']
        #if fw_spec['psuedo_potentials_path'] is not None:
        #    os.environ['PSP_PATH'] = fw_spec['psuedo_potentials_path']
        def Eng(abc, calcargs):
            #a,b,c = abc
            abc = np.array(abc)
            orig = dict_atoms(self['atoms'].copy())
            new = copy(orig)
            lattice_scale =  np.array([max(abs(a)) for a in atoms.cell])
            new.set_cell(np.multiply(abc/lattice_scale,new.cell.T).T,
                         scale_atoms=True)
            #print(new.posi)
            calc = SPARC(**calcargs)
            new.set_calculator(calc)
            E = new.get_potential_energy()
            del new, orig
            #MoO3.write('D-'+str(round(a,2))+'-'+str(round(b,2))+'-'+str(round(c,2))+'.traj')
            return E
        import numpy as np
        from scipy.optimize import fmin, minimize 
        from copy import copy
        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        #calc = SPARC(**parameter_dict)
        #atoms.set_calculator(calc)
        #calc.calculate()
        x0 = [max(abs(a)) for a in atoms.cell]
        #xopt = minimize(Eng,x0,args=(parameter_dict,))
        xopt = minimize(Eng, x0, args = (parameter_dict),method='Nelder-Mead',
                        options= {
                                    #'maxiter':150,
                                    'fatol':0.01,'xatol':0.01,
                                    #'adaptive': True
                                    }
                                        )
        x0 = [max(abs(a)) for a in atoms.cell]
        atoms.cell = np.multiply(xopt.x/x0,atoms.cell.T).T
        atoms.set_cell(np.multiply(xopt.x/x0,atoms.cell.T).T,
                       scale_atoms=True)
        calc = SPARC(**parameter_dict)
        atoms.set_calculator(calc)
        calc.calculate()

        if self['to_db'] == True:
            id_ = calc.calc_to_mongo(
                           )
            if self['identifier'] is not None:
                db = MongoDatabase(
                           ) 
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)

class OptimizeLatticeSPARC(Firework):
    def __init__(self, atoms,
                 parameters = default_parameters,
                 name = 'SPARC SCF',
                 sparc_command = None,
                 psuedo_potentials_path = None,
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 **kwargs
                 ):

        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(OptimizeLattice(atoms = atoms, parameter_dict = parameters,
             sparc_command = sparc_command,
             psuedo_potentials_path = psuedo_potentials_path,
             identifier = identifier,
             to_db = to_db))
        super(OptimizeLatticeSPARC, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

def get_sparc_lattice_optimizations(structures, parameters = default_parameters,
                                sparc_command = None,
                                to_db = True,
                                psuedo_potentials_path = None,
                                identifiers = None):
    fws = []
    if type(parameters) != list:  # If no list of parameters is given, use the same for all
        parameters = [parameters] * len(structures)
    if type(identifiers) != list and identifiers is not None:  # If no list of parameters is given, use the same for all
        identifiers = [identifiers] * len(structures)

    for struct, param, identifier in zip(structures, parameters,identifiers):
        name = struct.get_chemical_formula()
        fws.append(OptimizeLatticeSPARC(atoms_dict(struct),param,
                    sparc_command = sparc_command,
                    psuedo_potentials_path = psuedo_potentials_path,
                    identifier = identifier,
                    to_db = to_db))
    return Workflow(fws, name="{} tests wf, e.g.,".format(len(fws)))


@explicit_serialize
class OptimizeSparcASE(FiretaskBase):
    """
    Runs SPARC using ASE based on some input dictionary

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        
    """
    required_params = ['atoms', 'parameter_dict','to_db','identifier','db_file','fmax']
    optional_params = ['sparc_command']

    _fw_name = 'Run SPARC-X'
    def run_task(self,fw_spec):
        from ase.optimize import QuasiNewton
        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        calc = SPARC(**parameter_dict)
        atoms.set_calculator(calc)
        relax = QuasiNewton(atoms,logfile='opt.log',
                            trajectory='opt.traj',
                            restart='opt.pckl')
        relax.run(self['fmax'])
        if self['to_db'] == True:
            from json import load
            db_info = load(open(self['db_file'],'r'))
            id_ = calc.calc_to_mongo(**db_info )
            if self['identifier'] is not None:
                db = MongoDatabase(**db_info)
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)
        formula = atoms.get_chemical_formula()

        formula = atoms.get_chemical_formula()
        try:
            os.system('cp sprc-calc.Dens /gpfs/pace1/project/chbe-medford/medford-share/users/xlei38/sparc_w_print_executable/molecular_systems/density_files/' + formula + '.Dens')
        except:
            pass


class Sparc_Optimize_FW(Firework):
    def __init__(self, atoms,
                 parameters = default_parameters,
                 name = 'SPARC SCF',
                 sparc_command = None,
                 psuedo_potentials_path = None,
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 fmax = 0.02,
                 **kwargs
                 ):
        """
        Runs an SCF calculation in SPARC on the given structure

        Args:
            atoms (Atoms Object or Pymatgen Structure): Input structure
        
            name (str): Name for the Firework
        
            parameters (Dict): the input parameters for SPARC
        
            sparc_command (str): optional, a command used by ase to run SPARC.
                This can also be set with the $ASE_SPARC_COMMAND environment
                variable.

            psuedo_potentials_path (str): optional, a path to the pseudopotentials
                you'd like used. This can also be set using the $PSP_PATH environment
                Variable.

            db_file (str): Not implemented

            parents ([Firework]): Parents of this Firework
        """

        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(OptimizeSparcASE(atoms = atoms, parameter_dict = parameters,
             sparc_command = sparc_command,
             psuedo_potentials_path = psuedo_potentials_path,
             to_db = to_db,
             db_file = db_file,
             identifier = identifier,
             fmax = fmax))
        super(Sparc_Optimize_FW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

