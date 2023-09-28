import os
import numpy as np
from pathlib import Path
from ase.calculators.calculator import Calculator, FileIOCalculator, all_changes
import subprocess

from .io import SparcBundle
from .utils import _find_default_sparc, h2gpts, deprecated
from warnings import warn, warn_explicit
from .api import SparcAPI
import datetime
from ase.units import Bohr, Hartree, eV, GPa

# Below are a list of ASE-compatible calculator input parameters that are
# in Angstrom/eV units
# Ideas are taken from GPAW calculator
sparc_python_inputs = [
    "xc",
    "h",
    "kpts",
    "convergence",
    "gpts",
    "nbands",
]


defaultAPI = SparcAPI()


class SPARC(FileIOCalculator):
    # TODO: magmom should be a possible input
    implemented_properties = ["energy", "forces", "fermi", "stress"]
    name = "sparc"
    ase_objtype = "sparc_calculator"  # For JSON storage
    special_inputs = sparc_python_inputs

    # A "minimal" set of parameters that user can call plug-and-use
    # like atoms.calc = SPARC()
    # TODO: should we provide a minimal example for each system?
    default_params = {
        "xc": "pbe",
        "kpts": (1, 1, 1),
        "h": 0.25,  # Angstrom
    }

    def __init__(
        self,
        restart=None,
        directory=".",
        *,
        label=None,
        atoms=None,
        command=None,
        psp_dir=None,
        log="sparc.log",
        **kwargs,
    ):
        # Initialize the calculator but without restart.
        # Do not pass the label to the parent FileIOCalculator class to avoid issue
        # Handle old restart file separatedly since we rely on the sparc_bundle to work
        FileIOCalculator.__init__(
            self,
            restart=None,
            label=None,
            atoms=atoms,
            command=command,
            directory=directory,
            **kwargs,
        )

        # sparc bundle will set the label
        if label is None:
            label = "SPARC" if restart is None else None

        self.sparc_bundle = SparcBundle(
            directory=Path(self.directory),
            mode="w",
            atoms=self.atoms,
            label=label,
            psp_dir=psp_dir,
        )

        # Try restarting from an old calculation and set results
        self._restart(restart=restart)

        # Run a short test to return version of SPARC's binary
        # TODO: sparc_version should allow both read from results / short stdout
        self.sparc_version = self._detect_sparc_version()

        # Sanitize the kwargs by converting lower -- > upper
        # and perform type check
        # TODO: self.parameter should be the only entry
        self.valid_params, self.special_params = self._sanitize_kwargs(kwargs)
        self.log = self.directory / log if log is not None else None

    @property
    def directory(self):
        if hasattr(self, "sparc_bundle"):
            return Path(self.sparc_bundle.directory)
        else:
            return Path(self._directory)

    @directory.setter
    def directory(self, directory):
        if hasattr(self, "sparc_bundle"):
            self.sparc_bundle.directory = Path(directory)
        else:
            self._directory = Path(directory)
        return

    @property
    def label(self):
        """Rewrite the label from Calculator class, since we don't want to contain pathsep"""
        if hasattr(self, "sparc_bundle"):
            return self.sparc_bundle.label
        else:
            return getattr(self, "_label", None)

    @label.setter
    def label(self, label):
        """Rewrite the label from Calculator class,
        since we don't want to contain pathsep
        """
        label = str(label)
        if hasattr(self, "sparc_bundle"):
            self.sparc_bundle.label = self.sparc_bundle._make_label(label)
        else:
            self._label = label

    @property
    def sort(self):
        """Like Vasp calculator
        ASE atoms --> sort --> SPARC
        """
        return self.sparc_bundle.sorting["sort"]

    @property
    def resort(self):
        """Like Vasp calculator
        SPARC --> resort --> ASE atoms
        """
        return self.sparc_bundle.sorting["resort"]

    def check_state(self, atoms, tol=1e-9):
        """Updated check_state method.
        By default self.atoms (cached from output files) contains the initial_magmoms,
        so we add a zero magmoms to the atoms for comparison if it does not exist.

        reading a result from the .out file has only precision up to 10 digits
        """
        atoms_copy = atoms.copy()
        if "initial_magmoms" not in atoms_copy.arrays:
            atoms_copy.set_initial_magnetic_moments(
                [
                    0,
                ]
                * len(atoms_copy)
            )
        # First we check for default changes
        system_changes = FileIOCalculator.check_state(self, atoms_copy, tol=tol)
        # A few hard-written rules. Wrapping should only affect the position
        if "positions" in system_changes:
            atoms_copy.wrap()
            new_system_changes = FileIOCalculator.check_state(self, atoms_copy, tol=tol)
            # TODO: make sure such check only happens for PBC
            # the position is wrapped, accept as the same structure
            if "positions" not in new_system_changes:
                system_changes.remove("positions")
        return system_changes

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
                    "If this is not correct, "
                    "please manually set $ASE_SPARC_COMMAND or SPARC.command!"
                )
            self.command = command_env
        return f"{self.command} {extras}"

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        """Perform a calculation step"""
        # Check if the user accidentally provides atoms unit cell without vacuum

        if atoms and np.any(atoms.cell.cellpar()[:3] == 0):
            # TODO: choose a better error name
            msg = "Cannot setup SPARC calculation because at least one of the lattice dimension is zero!"
            if any([bc_ is False for bc_ in atoms.pbc]):
                msg += " Please add a vacuum in the non-periodic direction of your input structure."
            raise ValueError(msg)
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        self.execute()
        self.read_results()
        # Extra step, copy the atoms back to original atoms, if it's an
        # geopt or aimd calculation
        if ("geopt" in self.raw_results) or ("aimd" in self.raw_results):
            # Update the parent atoms
            atoms.set_positions(self.atoms.positions, apply_constraint=False)
            atoms.cell = self.atoms.cell
            atoms.constraints = self.atoms.constraints
            atoms.pbc = self.atoms.pbc
            # copy init magmom just to avoid check_state issue
            if "initial_magmoms" in self.atoms.arrays:
                atoms.set_initial_magnetic_moments(
                    self.atoms.get_initial_magnetic_moments()
                )

    def get_stress(self, atoms=None):
        """Warn user the dimensionality change when using stress"""
        if "stress_equiv" in self.results:
            raise NotImplementedError(
                "You're requesting stress in a low-dimensional system. Please use `calc.results['stress_equiv']` instead!"
            )
        return super().get_stress(atoms)

        # atoms = self.atoms.copy()

    # def update_atoms(self, atoms):
    #     """Update atoms after calculation if the positions are changed

    #     Idea taken from Vasp.update_atoms.
    #     """
    #     if (self.int_params['ibrion'] is not None
    #             and self.int_params['nsw'] is not None):
    #         if self.int_params['ibrion'] > -1 and self.int_params['nsw'] > 0:
    #             # Update atomic positions and unit cell with the ones read
    #             # from CONTCAR.
    #             atoms_sorted = read(self._indir('CONTCAR'))
    #             atoms.positions = atoms_sorted[self.resort].positions
    #             atoms.cell = atoms_sorted.cell

    # self.atoms = atoms  # Creates a copy

    def _check_input_exclusion(self, input_parameters, atoms=None):
        """Check if mutually exclusive parameters are provided

        The exclusion rules are taken from the SPARC manual and currently hard-coded.
        We may need to have a clever way to do the automatic rule conversion in API
        """
        # Rule 1: ECUT, MESH_SPACING, FD_GRID
        count = 0
        for key in ["ECUT", "MESH_SPACING", "FD_GRID"]:
            if key in input_parameters:
                count += 1
        if count > 1:
            # TODO: change to ExclusionParameterError
            raise ValueError(
                "ECUT, MESH_SPACING, FD_GRID cannot be specified simultaneously!"
            )

        # Rule 2: LATVEC_SCALE, CELL
        if ("LATVEC_SCALE" in input_parameters) and ("CELL" in input_parameters):
            # TODO: change to ExclusionParameterError
            raise ValueError(
                "LATVEC_SCALE and CELL cannot be specified simultaneously!"
            )

        # When the cell is provided via ase object, we will forbid user to provide
        # LATVEC, LATVEC_SCALE or CELL
        # TODO: make sure the rule makes sense for molecules
        if atoms is not None:
            if any([p in input_parameters for p in ["LATVEC", "LATVEC_SCALE", "CELL"]]):
                raise ValueError(
                    "When passing an ase atoms object, LATVEC, LATVEC_SCALE or CELL cannot be set simultaneously!"
                )

    def _check_minimal_input(self, input_parameters):
        """Check if the minimal input set is satisfied

        TODO: maybe we need to move the minimal set to class default
        """
        for param in ["EXCHANGE_CORRELATION", "KPOINT_GRID"]:
            if param not in input_parameters:
                # TODO: change to MissingParameterError
                raise ValueError(f"Parameter {param} is not provided.")
        # At least one from ECUT, MESH_SPACING and FD_GRID must be provided
        if not any(
            [param in input_parameters for param in ("ECUT", "MESH_SPACING", "FD_GRID")]
        ):
            raise ValueError(
                "You should provide at least one of ECUT, MESH_SPACING or FD_GRID."
            )

    def write_input(self, atoms, properties=[], system_changes=[]):
        """Create input files via SparcBundle"""
        # import pdb; pdb.set_trace()
        # print("Calling the properties: ", properties)
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        converted_params = self._convert_special_params(atoms=atoms)
        input_parameters = converted_params.copy()
        input_parameters.update(self.valid_params)

        # Make sure desired properties are always ensured, but we don't modify the user inputs
        if "forces" in properties:
            input_parameters["PRINT_FORCES"] = True

        if "stress" in properties:
            input_parameters["CALC_STRESS"] = True

        # TODO: detect if minimal values are set

        # TODO: system_changes ?

        # TODO: check parameter exclusion

        self._check_input_exclusion(input_parameters, atoms=atoms)
        self._check_minimal_input(input_parameters)

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
            input_parameters=input_parameters,
        )

    def execute(self):
        """Make the calculation. Note we probably need to use a better handling of background process!"""
        # TODO: add -socket?
        extras = f"-name {self.label}"
        command = self._make_command(extras=extras)
        self.print_sysinfo(command)

        # TODO: distinguish between normal process
        try:
            if self.log is not None:
                with open(self.log, "a") as fd:
                    self.proc = subprocess.run(
                        command, shell=True, cwd=self.directory, stdout=fd
                    )
            else:
                self.proc = subprocess.run(
                    command, shell=True, cwd=self.directory, stdout=None
                )
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        # We probably don't want to wait the
        errorcode = self.proc.returncode

        if errorcode > 0:
            msg = f"SPARC failed with command {command}" f"with error code {errorcode}"
            raise RuntimeError(msg)

        return

    @property
    def raw_results(self):
        return getattr(self.sparc_bundle, "raw_results", {})

    @raw_results.setter
    def raw_results(self, value):
        self.sparc_bundle.raw_results = value
        return

    def read_results(self):
        """Parse from the SparcBundle"""
        # TODO: try use cache?
        # self.sparc_bundle.read_raw_results()
        last = self.sparc_bundle.convert_to_ase(indices=-1, include_all_files=False)
        self.atoms = last.copy()
        self.results.update(last.calc.results)

        # self._extract_out_results()

    def _restart(self, restart=None):
        """Reload the input parameters and atoms from previous calculation.

        If self.parameters is already set, the parameters will not be loaded
        If self.atoms is already set, the atoms will be not be read
        """
        if restart is None:
            return
        reload_atoms = self.atoms is None
        reload_parameters = len(self.parameters) == 0

        self.read_results()
        if not reload_atoms:
            self.atoms = None
        if reload_parameters:
            self.parameters = self.raw_results["inpt"]["params"]

        if (not reload_parameters) or (not reload_atoms):
            warn(
                "Extra parameters or atoms are provided when restarting the SPARC calculator, "
                "previous results will be cleared."
            )
            self.results.clear()
            self.sparc_bundle.raw_results.clear()
        return

    def get_fermi_level(self):
        """Extra get-method for Fermi level, if calculated"""
        return self.results.get("fermi", None)

    def _detect_sparc_version(self):
        """Run a short sparc test to determine which sparc is used"""
        # TODO: complete the implementation
        return None

    def _sanitize_kwargs(self, kwargs):
        """Convert known parameters from"""
        # print(kwargs)
        # TODO: versioned validator
        validator = defaultAPI
        valid_params = {}
        special_params = self.default_params.copy()
        # TODO: how about overwriting the default parameters?
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
                    # TODO: should we raise exception instead?
                    warn(f"Input parameter {key} does not have a valid value!")
        return valid_params, special_params

    def _convert_special_params(self, atoms=None):
        """Convert ASE-compatible parameters to SPARC compatible ones
        parameters like `h`, `nbands` may need atoms information

        Special rules:
        h <--> gpts <--> FD_GRID, only when None of FD_GRID / ECUT or MESH_SPACING is provided
        """
        converted_sparc_params = {}
        validator = defaultAPI
        params = self.special_params.copy()

        # xc --> EXCHANGE_CORRELATION
        # TODO: more XC options
        if "xc" in params:
            xc = params.pop("xc")
            if xc.lower() == "pbe":
                converted_sparc_params["EXCHANGE_CORRELATION"] = "GGA_PBE"
            elif xc.lower() == "lda":
                converted_sparc_params["EXCHANGE_CORRELATION"] = "LDA_PW"

        # h --> gpts
        if "h" in params:
            if "gpts" in params:
                raise KeyError(
                    "h and gpts cannot be provided together in SPARC calculator!"
                )
            h = params.pop("h")
            if atoms is None:
                raise ValueError(
                    "Must have an active atoms object to convert h --> gpts!"
                )
            if any(
                [p in self.valid_params for p in ("FD_GRID", "ECUT", "MESH_SPACING")]
            ):
                warn(
                    "You have specified one of FD_GRID, ECUT or MESH_SPACING, "
                    "conversion of h to mesh grid is ignored."
                )
            else:
                # TODO: is there any limitation for parallelization?
                gpts = h2gpts(h, atoms.cell)
                params["gpts"] = gpts

        # gpts --> FD_GRID
        if "gpts" in params:
            gpts = params.pop("gpts")
            if validator.validate_input("FD_GRID", gpts):
                converted_sparc_params["FD_GRID"] = gpts
            else:
                # TODO: customize error
                raise ValueError(f"Input parameter gpts has invalid value {gpts}")

        # kpts
        if "kpts" in params:
            # TODO: how about accepting ASE's kpts setting?
            kpts = params.pop("kpts")
            if validator.validate_input("KPOINT_GRID", kpts):
                converted_sparc_params["KPOINT_GRID"] = kpts
            else:
                # TODO: customize error
                raise ValueError(f"Input parameter kpts has invalid value {kpts}")

        # nbands
        if "nbands" in params:
            # TODO: Check if the nbands are correct in current system
            # TODO: default $N_e/2 \\times 1.2 + 5$
            nbands = params.pop("nbands")
            if validator.validate_input("NSTATES", nbands):
                converted_sparc_params["NSTATES"] = nbands
            else:
                # TODO: customize error
                raise ValueError(f"Input parameter nbands has invalid value {nbands}")

        # convergence is a dict
        if "convergence" in params:
            convergence = params.pop("convergence")
            tol_e = convergence.get("energy", None)
            if tol_e:
                # TOL SCF: Ha / atom <--> energy tol: eV / atom
                converted_sparc_params["SCF_ENERGY_ACC"] = tol_e / Hartree

            # TODO: per AJ's suggestion, better change forces to relaxation
            tol_f = convergence.get("forces", None)
            if tol_f:
                # TOL SCF: Ha / Bohr <--> energy tol: Ha / Bohr
                converted_sparc_params["TOL_RELAX"] = tol_f / Hartree * Bohr

            tol_dens = convergence.get("density", None)
            if tol_dens:
                # TOL SCF: electrons / atom
                converted_sparc_params["TOL_PSEUDOCHARGE"] = tol_dens

            tol_stress = convergence.get("stress", None)
            if tol_stress:
                # TOL SCF: electrons / atom
                converted_sparc_params["TOL_RELAX_CELL"] = tol_stress / GPa

        return converted_sparc_params

    def print_sysinfo(self, command=None):
        """Record current runtime information"""
        now = datetime.datetime.now().isoformat()
        if command is None:
            command = self.command
        msg = (
            "\n" + "*" * 80 + "\n"
            f"SPARC program started by sparc-python-api at {now}\n"
            f"command: {command}\n"
        )
        if self.log is None:
            print(msg)
        else:
            with open(self.log, "a") as fd:
                print(msg, file=fd)

    ###############################################
    # Below are deprecated functions from v1
    ###############################################
    @deprecated("Please use SPARC.set instead for setting grid")
    def interpret_grid_input(self, atoms, **kwargs):
        return None

    @deprecated("Please use SPARC.set instead for setting kpoints")
    def interpret_kpoint_input(self, atoms, **kwargs):
        return None

    @deprecated("Please use SPARC.set instead for setting downsampling parameter")
    def interpret_downsampling_input(self, atoms, **kwargs):
        return None

    @deprecated("Please use SPARC.set instead for setting kpoint shift")
    def interpret_kpoint_shift(self, atoms, **kwargs):
        return None

    @deprecated("Please use SPARC.psp_dir instead")
    def get_pseudopotential_directory(self, pseudo_dir=None, **kwargs):
        return self.sparc_bundle.psp_dir

    # TODO: update method
    def get_nstates(self):
        return None

    @deprecated("Please set the variables separatedly")
    def setup_parallel_env(self):
        return None

    @deprecated("Please use SPARC._make_command instead")
    def generate_command(self):
        return self._make_command(f"-name {self.label}")

    # TODO: update the method!
    def estimate_memory(self, atoms=None, units="GB", **kwargs):
        """
        a function to estimate the amount of memory required to run
        the selected calculation. This function takes in **kwargs,
        but if none are passed in, it will fall back on the parameters
        input when the class was instantiated
        """
        conversion_dict = {
            "MB": 1e-6,
            "GB": 1e-9,
            "B": 1,
            "byte": 1,
            "KB": 1e-3,
        }
        if kwargs == {}:
            kwargs = self.parameters
        if atoms is None:
            atoms = self.atoms

        nstates = kwargs.get("NSTATES")
        if nstates is None:
            nstates = self.get_nstates(atoms=atoms, **kwargs)

        # some annoying code to figure out if it's a spin system
        spin_polarized = kwargs.get("nstates")
        if spin_polarized is not None:
            spin_polarized = int(spin_polarized)
        else:
            spin_polarized = 1
        if spin_polarized == 2:
            spin_factor = 2
        else:
            spin_factor = 1

        if "MESH_SPACING" in kwargs:
            kwargs["h"] = kwargs.pop("MESH_SPACING")
        npoints = np.product(self.interpret_grid_input(atoms, **kwargs))

        kpt_grid = self.interpret_kpoint_input(atoms, **kwargs)
        kpt_factor = np.ceil(np.product(kpt_grid) / 2)

        # this is a pretty generous over-estimate
        # TODO: check this function is working
        estimate = 5 * npoints * nstates * kpt_factor * spin_factor * 8  # bytes
        converted_estimate = estimate * conversion_dict[units]
        return converted_estimate

    # TODO: update method for static / geopt / aimd
    def get_scf_steps(self, include_uncompleted_last_step=False):
        raise NotImplemented

    @deprecated("Use SPARC.get_number_of_ionic_steps instead")
    def get_geometric_steps(self, include_uncompleted_last_step=False):
        raise NotImplemented

    def get_runtime(self):
        raise NotImplemented

    def get_fermi_level(self):
        raise NotImplemented

    # TODO: update method for static / geopt / aimd
    def get_scf_steps(self, include_uncompleted_last_step=False):
        raise NotImplemented

    @deprecated("Use SPARC.get_number_of_ionic_steps instead")
    def get_geometric_steps(self, include_uncompleted_last_step=False):
        raise NotImplemented

    def get_runtime(self):
        raise NotImplemented

    def get_fermi_level(self):
        raise NotImplemented

    @deprecated
    def concatinate_output(self):
        raise DeprecationWarning("Functionality moved in sparc.SparcBundle.")

    @deprecated
    def read_line(self, **kwargs):
        raise DeprecationWarning(
            "Parsers for individual files have been moved to sparc.sparc_parsers module"
        )

    @deprecated
    def parse_output(self, **kwargs):
        raise DeprecationWarning("Use SPARC.read_results for parsing results!")

    @deprecated
    def parse_relax(self, *args, **kwargs):
        raise DeprecationWarning("Use SPARC.read_results for parsing results!")

    @deprecated
    def parse_MD(self, *args, **kwargs):
        raise DeprecationWarning("Use SPARC.read_results for parsing results!")

    @deprecated
    def parse_input_args(self, input_block):
        raise DeprecationWarning("Use SPARC.set for argument handling!")

    @deprecated
    def recover_index_order_from_ion_file(self, label):
        raise DeprecationWarning(
            "Use SPARC.sort and SPARC.resort for atomic index sorting!"
        )

    @deprecated
    def atoms_dict(self, *args, **kwargs):
        raise DeprecationWarning("")

    @deprecated
    def dict_atoms(self, *args, **kwargs):
        raise DeprecationWarning("")
