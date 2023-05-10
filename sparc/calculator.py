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
import subprocess

from .sparc_io_bundle import SparcBundle
from .utils import _find_default_sparc
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
    name = "sparc"
    ase_objtype = "sparc_calculator"  # For JSON storage
    special_inputs = sparc_python_inputs

    def __init__(
        self,
        restart=None,
        directory=".",
        *,
        label=None,
        atoms=None,
        command=None,
        psp_dir=None,
        **kwargs,
    ):
        FileIOCalculator.__init__(
            self, restart=restart, label=label, atoms=atoms, command=command, **kwargs
        )

        self.directory = Path(directory)

        # SparcBundle.__init__(        #                        directory=directory,
        #                        mode="w",
        #                        atoms=atoms,
        #                        label=label,
        #                        psp_dir=psp_dir)

        # TODO: change label?
        self.label = label if label is not None else "SPARC"
        self.sparc_bundle = SparcBundle(
            directory=self.directory,
            mode="w",
            atoms=atoms,
            label=self.label,
            psp_dir=psp_dir,
        )
        # Run a short test to return version of SPARC's binary
        self.sparc_version = self._detect_sparc_version()
        # Sanitize the kwargs by converting lower -- > upper
        # and perform type check
        self.valid_params, self.special_params = self._sanitize_kwargs(kwargs)

    def _make_command(self, extras=""):
        """Use $ASE_SPARC_COMMAND or self.command to determine the command
        as a last resort, if `sparc` exists in the PATH, use that information

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
                sparc_exe, mpi_exe, num_cores = _find_default_sparc()
                if sparc_exe is None:
                    raise EnvironmentError(
                        "Cannot find your sparc setup via $ASE_SPARC_COMMAND, SPARC.command, or "
                        "infer from your $PATH. Please refer to the manual!"
                    )
                if mpi_exe is not None:
                    command_env = f"{mpi_exe} -n {num_cores} {sparc_exe}"
                else:
                    command_env = f"{sparc_exe}"
                warn(
                    f"Your sparc command is inferred to be {command_env}, "
                    "If this is not correct, please manually set $ASE_SPARC_COMMAND or SPARC.command!"
                )
            self.command = command_env
        return f"{self.command} {extras}"

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        self.execute()
        self.read_results()

    def write_input(self, atoms, properties=None, system_changes=None):
        """Create input files via SparcBundle"""
        # import pdb; pdb.set_trace()
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        self.sparc_bundle._write_ion_and_inpt(
            atoms=atoms,
            label=self.label,
            # Pass the rest parameters from calculator!
            direct=False,
            sort=True,
            ignore_constraints=False,
            wrap=False,
            # Below are the parameters from v1
            # scaled -> direct, ignore_constraints --> not add_constraints
            scaled=False,
            add_constraints=True,
            copy_psp=True,
            comment="",
            input_parameters=self.valid_params,
        )

    def execute(self):
        """Make the calculation. Note we probably need to use a better handling of background process!"""
        # TODO: add -socket?
        extras = f"-name {self.label}"
        command = self._make_command(extras=extras)

        # TODO: distinguish between normal process
        try:
            self.proc = subprocess.run(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        # We probably don't want to wait the
        errorcode = self.proc.returncode

        if errorcode > 0:
            msg = f"SPARC failed with command {command}" f"with error code {errorcode}"
            raise CalculationFailed(msg)

        return

    def read_results(self):
        """Parse from the SparcBundle"""
        pass

    def _detect_sparc_version(self):
        """Run a short sparc test to determine which sparc is used"""
        # TODO: complete the implementation
        return None

    def _sanitize_kwargs(self, kwargs):
        """Convert known parameters from"""
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
                    warn(
                        f"Parameter {key} (case-insentivie) appears multiple times in the calculator setup!"
                    )
                if validator.validate_input(key, value):
                    valid_params[key] = value
                else:
                    # TODO: helper information
                    warn(f"Input parameter {key} does not have a valid value!")
        return valid_params, special_params
