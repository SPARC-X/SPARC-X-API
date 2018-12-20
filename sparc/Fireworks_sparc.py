from .sparc import SPARC, default_parameters
from .utilities import atoms_dict, dict_atoms
import os
from fireworks import explicit_serialize, FiretaskBase, FWAction, Firework, Workflow
from ase.atoms import Atoms
from ase.atom import Atom
from .mongo import mongo_atoms_doc, mongo_doc_atoms
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
    required_params = ['atoms', 'parameter_dict']
    optional_params = ['sparc_command']

    _fw_name = 'Write Sparc Input From Dict'
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
        calc.calc_to_mongo( 
                           host='ds153593.mlab.com',
                           port=53593,
                           database='comer-test-db',
                           user='AJ',
                           password='a123456',)

class Sparc_SCF_FW(Firework):
    def __init__(self, atoms, 
                 parameters = default_parameters,  
                 name = 'SPARC SCF',
                 sparc_command = None,
                 psuedo_potentials_path = None,
                 db_file = None,
                 parents = None,
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
             psuedo_potentials_path = psuedo_potentials_path))
        super(Sparc_SCF_FW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

def get_sparc_convergence_tests(structures, parameters = default_parameters,
                                sparc_command = None,
                                psuedo_potentials_path = None):
    fws = []
    if type(parameters) != list:  # If no list of parameters is given, use the same for all
        parameters = [parameters] * len(structures) 
    for struct, param in zip(structures, parameters):
        name = struct.get_chemical_formula()
        fws.append(Sparc_SCF_FW(atoms_dict(struct),param,
                    sparc_command = sparc_command,
                    psuedo_potentials_path = psuedo_potentials_path))
    return Workflow(fws, name="{} surfaces wf, e.g., {}".format(len(fws), fws[0].name))
