import numpy as np
import os
from pathlib import Path
from ase.calculators.calculator import FileIOCalculator, all_changes
from ase.units import Bohr, Hartree, fs, GPa
from ase.calculators.calculator import CalculatorError, CalculatorSetupError
from ase.calculators.calculator import CalculationFailed, SCFError, ReadError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import compare_atoms
from ase.atoms import Atoms

from .sparc_io_bundle import SparcBundle
from warnings import warn

# Below are a list of ASE-compatible calculator input parameters that are
# in Angstrom/eV units
# Ideas are taken from GPAW calculator
sparc_python_inputs = ["xc", "h", "kpts", "convergence", "gpts", "nbands", "encut"]

from .inputs import SparcInputs

defaultAPI = SparcInputs()


class SPARC(SparcBundle, FileIOCalculator):
    # TODO: magmom should be a possible input
    implemented_properties = ["energy", "forces", "fermi", "stress"]
    name = 'sparc'
    ase_objtype = 'sparc_calculator'  # For JSON storage
    special_inputs = sparc_python_inputs

    def __init__(self,
                 restart=None,
                 directory=".",
                 *,
                 label=None,
                 atoms=None,
                 command=None,
                 psp_dir=None,
                 **kwargs):
        
        FileIOCalculator.__init__(self,
                                  restart=restart,
                                  label=label,
                                  atoms=atoms,
                                  command=command,
                                  **kwargs)

        self.directory = Path(self.directory)

        SparcBundle.__init__(        #                        directory=directory,
        #                        mode="w",
        #                        atoms=atoms,
        #                        label=label,
        #                        psp_dir=psp_dir)
        # self.sparc_bundle = SparcBundle(
        #                        directory=directory,
        #                        mode="w",
        #                        atoms=atoms,
        #                        label=label,
        #                        psp_dir=psp_dir)
        # Run a short test to return version of SPARC's binary
        self.sparc_version = self._detect_sparc_version()
        # Sanitize the kwargs by converting lower -- > upper
        # and perform type check
        self.valid_params, self.special_params = self._sanitize_kwargs(kwargs)
        


    def _make_command(self, extras=""):
        """Use $ASE_SPARC_COMMAND or self.command to determine the command

        Extras will add additional arguments to the self.command,
        e.g. -name, -socket etc
        """
        if isinstance(extras, (list, tuple)):
            extras = " ".join(extras)
        else:
            extras = extras.strip()
        if self.command is None:
            command_env = os.environ.get("ASE_SPARC_COMMAND", None)
            if command_env is None:
                raise ValueError("Please set either $ASE_SPARC_COMMAND or SPARC.command!")
            self.command = command_env
        return self.command + extras
    

    def _detect_sparc_version(self):
        """Run a short sparc test to determine which sparc is used"""
        #TODO: complete the implementation
        return None
        

    def _sanitize_kwargs(self, kwargs):
        """Convert known parameters from
        """
        print(kwargs)
        # TODO: versioned validator
        validator = defaultAPI
        valid_params = {}
        special_params = {}
        # SPARC API is case insensitive
        for key, value in kwargs.items():
            if key in self.special_inputs:
                special_params[key] = value
            else:
                key = key.upper()
                if key in valid_params:
                    warn(f"Parameter {key} (case-insentivie) appears multiple times in the calculator setup!")
                if validator.validate_input(key, value):
                    valid_params[key] = value
                else:
                    # TODO: helper information
                    warn(f"Input parameter {key} does not have a valid value!")
        return valid_params, special_params

    def write_input(self)
        
