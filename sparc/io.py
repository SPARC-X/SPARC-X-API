"""Providing a new bundled SPARC file format
"""
import os
import re
from pathlib import Path
from warnings import warn

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator

# various io formatters
from .api import SparcAPI
from .common import psp_dir as default_psp_dir
from .download_data import is_psp_download_complete
from .sparc_parsers.aimd import _read_aimd
from .sparc_parsers.atoms import atoms_to_dict, dict_to_atoms
from .sparc_parsers.geopt import _read_geopt
from .sparc_parsers.inpt import _read_inpt, _write_inpt
from .sparc_parsers.ion import _read_ion, _write_ion
from .sparc_parsers.out import _read_out
from .sparc_parsers.pseudopotential import copy_psp_file, parse_psp8_header
from .sparc_parsers.static import _add_cell_info, _read_static
from .utils import deprecated, locate_api, string2index

# from .sparc_parsers.ion import read_ion, write_ion
defaultAPI = locate_api()


class SparcBundle:
    """Provide access to a calculation folder of SPARC as a simple bundle

    The bundle can be optionally named as .sparc following the ASE's
    .bundle format

    Currently the write method only supports 1 image, while read method support reading
    atoms results in following conditions

    1) No calculation (minimal): .ion + .inpt file --> 1 image
    2) Single point calculation: .ion + .inpt + .out + .static --> 1
    image with calc
    3) Multiple SP calculations: chain all
    .out{digits} and .static{digitis} outputs 4) Relaxation: read from
    .geopt and .out (supporting chaining) 5) AIMD: read from .aimd and
    .out (support chaining)


    Attributes:
        directory (Path): Path to the directory containing SPARC files.
        mode (str): File access mode ('r', 'w', or 'a').
        label (str): Name of the main SPARC file.
        init_atoms (Atoms): Initial atomic configuration.
        init_inputs (dict): Initial input parameters.
        psp_data (dict): Pseudopotential data.
        raw_results (dict): Raw results from SPARC calculations.
        psp_dir (Path): Directory containing pseudopotentials.
        sorting (list): Sort order for atoms.
        last_image (int): Index of the last image in a series of calculations.
        validator (SparcAPI): API validator for SPARC calculations.

    Methods:
        __find_psp_dir(psp_dir=None): Finds the directory for SPARC pseudopotentials.
        _find_files(): Finds all files matching the bundle label.
        _make_label(label=None): Infers or sets the label for the SPARC bundle.
        _indir(ext, label=None, occur=0, d_format="{:02d}"): Finds a file with a specific extension in the bundle.
        _read_ion_and_inpt(): Reads .ion and .inpt files together.
        _write_ion_and_inpt(): Writes .ion and .inpt files to the bundle.
        _read_results_from_index(index, d_format="{:02d}"): Reads results from a specific calculation index.
        _make_singlepoint(calc_results, images, raw_results): Converts results and images to SinglePointDFTCalculators.
        _extract_static_results(raw_results, index=":"): Extracts results from static calculations.
        _extract_geopt_results(raw_results, index=":"): Extracts results from geometric optimization calculations.
        _extract_aimd_results(raw_results, index=":"): Extracts results from AIMD calculations.
        convert_to_ase(index=-1, include_all_files=False, **kwargs): Converts raw results to ASE Atoms with calculators.
        read_raw_results(include_all_files=False): Parses all files in the bundle and merges results.
        read_psp_info(): Parses pseudopotential information from the inpt file.
    """

    psp_env = ["SPARC_PSP_PATH", "SPARC_PP_PATH"]

    def __init__(
        self,
        directory,
        mode="r",
        atoms=None,
        label=None,
        psp_dir=None,
        validator=defaultAPI,
    ):
        """
        Initializes a SparcBundle for accessing SPARC calculation data.

        Args:
            directory (str or Path): The path to the directory containing the SPARC files.
            mode (str, optional): The file access mode. Can be 'r' (read), 'w' (write), or 'a' (append). Defaults to 'r'.
            atoms (Atoms, optional): The initial atomic configuration. Only relevant in write mode.
            label (str, optional): A custom label for the bundle. If None, the label is inferred from the directory or files.
            psp_dir (str or Path, optional): Path to the directory containing pseudopotentials. If None, the path is inferred.
            validator (SparcAPI, optional): An instance of SparcAPI for validating and parsing SPARC parameters. Defaults to a default SparcAPI instance.

        Raises:
            AssertionError: If an invalid mode is provided.
            ValueError: If multiple .ion files are found and no label is specified.
            Warning: If no .ion file is found in read-mode, or illegal characters are in the label.
        """
        self.directory = Path(directory)
        self.mode = mode.lower()
        assert self.mode in (
            "r",
            "w",
            "a",
        ), f"Invalid mode {self.mode}! Must one of 'r', 'w' or 'a'"
        self.label = self._make_label(label)
        self.init_atoms = atoms.copy() if atoms is not None else None
        self.init_inputs = {}
        self.psp_data = {}
        self.raw_results = {}
        self.psp_dir = self.__find_psp_dir(psp_dir)
        # Sorting should be consistent across the whole bundle!
        self.sorting = None
        self.last_image = -1
        self.validator = validator

    def _find_files(self):
        """Find all files matching '{label}.*'"""
        return list(self.directory.glob(f"{self.label}.*"))

    def _make_label(self, label=None):
        """Infer the label from the bundle

        Special cases if label is None:
        1. read mode --> get the ion file name
        2. write mode --> infer from the directory

        Arguments:
            label (str or None): Label to be used to write the .ion, .inpt files
        """
        prefix = self.directory.resolve().with_suffix("").name

        illegal_chars = '\\/:*?"<>|'
        if label is not None:
            label_ = label
        elif self.mode == "w":
            label_ = prefix
        else:
            # read
            match_ion = list(self.directory.glob("*.ion"))
            if len(match_ion) > 1:
                raise ValueError(
                    "Cannot read sparc bundle with multiple ion files without specifying the label!"
                )
            elif len(match_ion) == 1:
                label_ = match_ion[0].name.split(".")[0]
            else:
                # No file found, possibly an empty bundle
                warn("No .ion file found in the read-mode bundle.")
                label_ = prefix

        if any([c in label_ for c in illegal_chars]):
            warn(
                f"Label name {label_} contains illegal characters! I'll make it 'SPARC'"
            )
            label_ = "SPARC"
        return label_

    def __find_psp_dir(self, psp_dir=None):
        """Use environmental variable to find the directory for SPARC
        pseudopotentials

        Searching priority:
        1. User defined psp_dir
        2. $SPARC_PSP_PATH
        3. $SPARC_PP_PATH
        4. psp bundled with sparc-api

        Arguments:
            psp_dir (str or PosixPath or None): the specific directory to search the psp files.
                                                Each element can only have 1 psp file under psp_dir
        Returns:
            PosixPath: Location of psp files
        """
        if psp_dir is not None:
            return Path(psp_dir)
        else:
            for var in self.psp_env:
                env_psp_dir = os.environ.get(var, None)
                if env_psp_dir:
                    return Path(env_psp_dir)
            # At this point, we try to use the psp files bundled with sparc
            if is_psp_download_complete(default_psp_dir):
                return default_psp_dir
            else:
                warn(
                    (
                        "PSP directory bundled with SPARC-X-API is broken! "
                        "Please use `sparc.download_data` to re-download them!"
                    )
                )

            # Not found
            if self.mode == "w":
                warn(
                    (
                        "No pseudopotential searching path was set and "
                        "neither of $SPARC_PSP_PATH nor $SPARC_PP_PATH is set.\n"
                        "Please explicitly provide the pseudopotentials parameter when writing the sparc bundle."
                    )
                )
            return None

    def _indir(self, ext, label=None, occur=0, d_format="{:02d}"):
        """Find the file with {label}.{ext} under current dir,
        if label is None, use the default

        Arguments:
            ext (str): Extension of file, e.g. '.ion' or 'ion'
            label (str or None): Label for the file. If None, use the parent directory name for searching
            occur (int): Occurance index of the file, if occur > 0, search for files with suffix like 'SPARC.out_01'
            d_format (str): Format for the index

        Returns:
            PosixPath: Path to the target file under self.directory
        """
        label = self.label if label is None else label
        if not ext.startswith("."):
            ext = "." + ext
        if occur == 0:
            target = self.directory / f"{label}{ext}"
        else:
            target = self.directory / f"{label}{ext}_{d_format.format(occur)}"
        return target

    def _read_ion_and_inpt(self):
        """Read the ion and inpt files together to obtain basic atomstic data.

        Returns:
            Atoms: atoms object from .ion and .inpt file
        """
        f_ion, f_inpt = self._indir(".ion"), self._indir(".inpt")
        ion_data = _read_ion(f_ion, validator=self.validator)
        inpt_data = _read_inpt(f_inpt, validator=self.validator)
        merged_data = {**ion_data, **inpt_data}
        return dict_to_atoms(merged_data)

    def _write_ion_and_inpt(
        self,
        atoms=None,
        label=None,
        direct=False,
        sort=True,
        ignore_constraints=False,
        wrap=False,
        # Below are the parameters from v1
        # scaled -> direct, ignore_constraints --> not add_constraints
        scaled=False,
        add_constraints=True,
        copy_psp=False,
        comment="",
        input_parameters={},
        # Parameters that do not require type conversion
        **kwargs,
    ):
        """Write the ion and inpt files to a bundle. This method only
        supports writing 1 image.  If input_parameters are empty,
        there will only be .ion writing the positions and .inpt
        writing a minimal cell information

        Args:
            atoms (Atoms, optional): The Atoms object to write. If None, uses initialized atoms associated with SparcBundle.
            label (str, optional): Custom label for the written files.
            direct (bool, optional): If True, writes positions in direct coordinates.
            sort (bool, optional): If True, sorts atoms before writing.
            ignore_constraints (bool, optional): If True, ignores constraints on atoms.
            wrap (bool, optional): If True, wraps atoms into the unit cell.
            **kwargs: Additional keyword arguments for writing.

        Raises:
            ValueError: If the bundle is not in write mode.
        """
        if self.mode != "w":
            raise ValueError(
                "Cannot write input files while sparc bundle is opened in read or append mode!"
            )
        os.makedirs(self.directory, exist_ok=True)
        atoms = self.atoms.copy() if atoms is None else atoms.copy()
        pseudopotentials = kwargs.pop("pseudopotentials", {})

        if sort:
            if self.sorting is not None:
                old_sort = self.sorting.get("sort", None)
                if old_sort:
                    sort = old_sort

        data_dict = atoms_to_dict(
            atoms,
            direct=direct,
            sort=sort,
            ignore_constraints=ignore_constraints,
            psp_dir=self.psp_dir,
            pseudopotentials=pseudopotentials,
        )
        merged_inputs = input_parameters.copy()
        merged_inputs.update(kwargs)
        data_dict["inpt"]["params"].update(merged_inputs)

        # If copy_psp, change the PSEUDO_POT field and copy the files
        if copy_psp:
            for block in data_dict["ion"]["atom_blocks"]:
                if "PSEUDO_POT" in block:
                    origin_psp = block["PSEUDO_POT"]
                    target_dir = self.directory
                    target_fname = copy_psp_file(origin_psp, target_dir)
                    block["PSEUDO_POT"] = target_fname

        _write_ion(self._indir(".ion"), data_dict, validator=self.validator)
        _write_inpt(self._indir(".inpt"), data_dict, validator=self.validator)
        # Update the sorting information
        ion_dict = _read_ion(self._indir(".ion"))["ion"]
        self.sorting = ion_dict.get("sorting", None)
        return

    def read_raw_results(self, include_all_files=False):
        """Parse all files using the given self.label.
        The results are merged dict from all file formats

        Arguments:
            include_all_files (bool): Whether to include output files with different suffices
                                      If true: include all files (e.g. SPARC.out, SPARC.out_01,
                                      SPARC.out_02, etc).
        Returns:
            dict or List: Dict containing all raw results. Only some of them will appear in the calculator's results

        Sets:
            self.raw_results (dict or List): the same as the return value

        #TODO: @TT 2024-11-01 allow accepting indices
        #TODO: @TT last_image is a bad name, it should refer to the occurance of images
               the same goes with num_calculations
        """
        # Find the max output index
        out_files = self.directory.glob(f"{self.label}.out*")
        valid_out_files = [
            f
            for f in out_files
            if (re.fullmatch(r"^\.out(?:_\d+)?$", f.suffix) is not None)
        ]
        # Combine and sort the file lists
        last_out = sorted(valid_out_files, reverse=True)
        # No output file, only ion / inpt
        if len(last_out) == 0:
            self.last_image = -1
        else:
            suffix = last_out[0].suffix
            if suffix == ".out":
                self.last_image = 0
            else:
                self.last_image = int(suffix.split("_")[1])
        self.num_calculations = self.last_image + 1

        # Always make sure ion / inpt results are parsed regardless of actual calculations
        if include_all_files:
            if self.num_calculations > 0:
                results = [
                    self._read_results_from_index(index)
                    for index in range(self.num_calculations)
                ]
            else:
                results = [self._read_results_from_index(self.last_image)]
        else:
            results = self._read_results_from_index(self.last_image)

        self.raw_results = results

        if include_all_files:
            init_raw_results = self.raw_results[0]
        else:
            init_raw_results = self.raw_results.copy()

        self.init_atoms = dict_to_atoms(init_raw_results)
        self.init_inputs = {
            "ion": init_raw_results["ion"],
            "inpt": init_raw_results["inpt"],
        }
        self.psp_data = self.read_psp_info()
        return self.raw_results

    def _read_results_from_index(self, index, d_format="{:02d}"):
        """Read the results from one calculation index, and return a
        single raw result dict, e.g. for index=0 --> .static
        and index=1 --> .static_01.

        Arguments:
            index (int): Index of image to return the results
            d_format (str): Format for the index suffix

        Returns:
            dict: Results for single image

        #TODO: @TT should we call index --> occurance?

        """
        results_dict = {}

        for ext in ("ion", "inpt"):
            f = self._indir(ext, occur=0)
            if f.is_file():
                data_dict = globals()[f"_read_{ext}"](f)
                results_dict.update(data_dict)
        for ext in ("geopt", "static", "aimd", "out"):
            f = self._indir(ext, occur=index, d_format=d_format)
            if f.is_file():
                data_dict = globals()[f"_read_{ext}"](f)
                results_dict.update(data_dict)

        # Must have files: ion, inpt
        if ("ion" not in results_dict) or ("inpt" not in results_dict):
            raise RuntimeError(
                "Either ion or inpt files are missing from the bundle! "
                "Your SPARC calculation may be corrupted."
            )

        # Copy the sorting information, if not existing
        sorting = results_dict["ion"].get("sorting", None)
        if sorting is not None:
            if self.sorting is None:
                self.sorting = sorting
            else:
                # Compare stored sorting
                assert (tuple(self.sorting["sort"]) == tuple(sorting["sort"])) and (
                    tuple(self.sorting["resort"]) == tuple(sorting["resort"])
                ), "Sorting information changed!"
        return results_dict

    def convert_to_ase(self, index=-1, include_all_files=False, **kwargs):
        """Read the raw results from the bundle and create atoms with
        single point calculators

        Arguments:
            index (int or str): Index or slice of the image(s) to convert. Uses the same format as ase.io.read
            include_all_files (bool): If true, also read results with indexed suffices

        Returns:
            Atoms or List[Atoms]: ASE-atoms or images with single point results

        """
        # Convert to images!
        # TODO: @TT 2024-11-01 read_raw_results should implement a more
        # robust behavior handling index, as it is the entry point for all
        rs = self.read_raw_results(include_all_files=include_all_files)
        if isinstance(rs, dict):
            raw_results = [rs]
        else:
            raw_results = list(rs)
        res_images = []
        for entry in raw_results:
            if "static" in entry:
                calc_results, images = self._extract_static_results(entry, index=":")
            elif "geopt" in entry:
                calc_results, images = self._extract_geopt_results(entry, index=":")
            elif "aimd" in entry:
                calc_results, images = self._extract_aimd_results(entry, index=":")
            else:
                calc_results, images = None, [self.init_atoms.copy()]

            if images is not None:
                if calc_results is not None:
                    images = self._make_singlepoint(calc_results, images, entry)
                res_images.extend(images)

        if isinstance(index, int):
            return res_images[index]
        else:
            return res_images[string2index(index)]

    def _make_singlepoint(self, calc_results, images, raw_results):
        """Convert a calculator dict and images of Atoms to list of
        SinglePointDFTCalculators

        The calculator also takes parameters from ion, inpt that exist
        in self.raw_results.

        Arguments:
            calc_results (List): Calculation results for all images
            images (List): Corresponding Atoms images
            raw_results (List): Full raw results dict to obtain additional information

        Returns:
            List(Atoms): ASE-atoms images with single point calculators attached

        """
        converted_images = []
        for res, _atoms in zip(calc_results, images):
            atoms = _atoms.copy()
            sp = SinglePointDFTCalculator(atoms)
            # Res can be empty at this point, leading to incomplete calc
            sp.results.update(res)
            sp.name = "sparc"
            sp.kpts = raw_results["inpt"].get("params", {}).get("KPOINT_GRID", None)
            # There may be a better way handling the parameters...
            sp.parameters = raw_results["inpt"].get("params", {})
            sp.raw_parameters = {
                "ion": raw_results["ion"],
                "inpt": raw_results["inpt"],
            }
            atoms.calc = sp
            converted_images.append(atoms)
        return converted_images

    def _extract_static_results(self, raw_results, index=":"):
        """Extract the static calculation results and atomic
        structure(s) Returns: calc_results: dict with at least energy
        value atoms: ASE atoms object The priority is to parse
        position from static file first, then fallback from ion + inpt

        Note: make all energy / forces resorted!

        Arguments:
            raw_results (dict): Raw results parsed from self.read_raw_results
            index (str or int): Index or slice of images

        Returns:
            List[results], List[Atoms]

        """
        static_results = raw_results.get("static", [])
        calc_results = []
        # Use extra lattice information to construct the positions
        cell = self.init_atoms.cell
        # import pdb; pdb.set_trace()
        static_results = _add_cell_info(static_results, cell)

        if isinstance(index, int):
            _images = [static_results[index]]
        elif isinstance(index, str):
            _images = static_results[string2index(index)]

        ase_images = []
        for static_results in _images:
            partial_results = {}
            if "free energy" in static_results:
                partial_results["energy"] = static_results["free energy"]
                partial_results["free energy"] = static_results["free energy"]

            if "forces" in static_results:
                partial_results["forces"] = static_results["forces"][self.resort]

            if "atomic_magnetization" in static_results:
                partial_results["magmoms"] = static_results["atomic_magnetization"][
                    self.resort
                ]

            if "net_magnetization" in static_results:
                partial_results["magmom"] = static_results["net_magnetization"]

            if "stress" in static_results:
                partial_results["stress"] = static_results["stress"]

            if "stress_equiv" in static_results:
                partial_results["stress_equiv"] = static_results["stress_equiv"]

            atoms = self.init_atoms.copy()
            # import pdb; pdb.set_trace()
            if "atoms" in static_results:
                atoms_dict = static_results["atoms"]

                # The socket mode case. Reset all cell and positions
                # Be careful,
                if "lattice" in static_results:
                    lat = static_results["lattice"]
                    atoms.set_cell(lat, scale_atoms=False)
                    if "coord" not in atoms_dict:
                        raise KeyError(
                            "Coordination conversion failed in socket static output!"
                        )
                    atoms.set_positions(
                        atoms_dict["coord"][self.resort], apply_constraint=False
                    )
                else:  # Do not change cell information (normal static file)
                    if "coord_frac" in atoms_dict:
                        atoms.set_scaled_positions(
                            atoms_dict["coord_frac"][self.resort]
                        )
                    elif "coord" in atoms_dict:
                        atoms.set_positions(
                            atoms_dict["coord"][self.resort], apply_constraint=False
                        )
            ase_images.append(atoms)
            calc_results.append(partial_results)
        return calc_results, ase_images

    def _extract_geopt_results(self, raw_results, index=":"):
        """Extract the static calculation results and atomic
        structure(s) Returns: calc_results: dict with at least energy
        value atoms: ASE atoms object The priority is to parse
        position from static file first, then fallback from ion + inpt

        Arguments:
            raw_results (dict): Raw results parsed from self.read_raw_results
            index (str or int): Index or slice of images

        Returns:
            List[results], List[Atoms]

        """
        # print("RAW_RES:  ", raw_results)
        geopt_results = raw_results.get("geopt", [])
        calc_results = []
        if len(geopt_results) == 0:
            warn(
                "Geopt file is empty! This is not an error if the calculation is continued from restart. "
            )
            return None, None

        if isinstance(index, int):
            _images = [geopt_results[index]]
        elif isinstance(index, str):
            _images = geopt_results[string2index(index)]

        ase_images = []
        for result in _images:
            atoms = self.init_atoms.copy()
            partial_result = {}
            if "energy" in result:
                partial_result["energy"] = result["energy"]
                partial_result["free energy"] = result["energy"]

            if "forces" in result:
                partial_result["forces"] = result["forces"][self.resort]

            if "stress" in result:
                partial_result["stress"] = result["stress"]

            # Modify the atoms copy
            if "positions" in result:
                atoms.set_positions(
                    result["positions"][self.resort], apply_constraint=False
                )
                if "ase_cell" in result:
                    atoms.set_cell(result["ase_cell"])
            else:
                # For geopt and RELAX=2 (cell relaxation),
                # the positions may not be written in .geopt file
                relax_flag = raw_results["inpt"]["params"].get("RELAX_FLAG", 0)
                if relax_flag != 2:
                    raise ValueError(
                        ".geopt file missing positions while RELAX!=2. "
                        "Please check your setup ad output files."
                    )
                if "ase_cell" not in result:
                    raise ValueError(
                        "Cannot recover positions from .geopt file due to missing cell information. "
                        "Please check your setup ad output files."
                    )
                atoms.set_cell(result["ase_cell"], scale_atoms=True)

            # Unlike low-dimensional stress in static calculations, we need to convert
            # stress_1d stress_2d to stress_equiv using the non-period cell dimension(s)
            # This has to be done when the actual cell information is loaded
            if "stress_1d" in result:
                stress_1d = result["stress_1d"]
                assert (
                    np.count_nonzero(atoms.pbc) == 1
                ), "Dimension of stress and PBC mismatch!"
                for i, bc in enumerate(atoms.pbc):
                    if not bc:
                        stress_1d /= atoms.cell.cellpar()[i]
                stress_equiv = stress_1d
                partial_result["stress_equiv"] = stress_equiv

            if "stress_2d" in result:
                stress_2d = result["stress_2d"]
                assert (
                    np.count_nonzero(atoms.pbc) == 2
                ), "Dimension of stress and PBC mismatch!"
                for i, bc in enumerate(atoms.pbc):
                    if not bc:
                        stress_2d /= atoms.cell.cellpar()[i]
                stress_equiv = stress_2d
                partial_result["stress_equiv"] = stress_equiv

            calc_results.append(partial_result)
            ase_images.append(atoms)

        return calc_results, ase_images

    def _extract_aimd_results(self, raw_results, index=":"):
        """Extract energy / forces from aimd results

        For calculator, we only need the last image

        We probably want more information for the AIMD calculations,
        but I'll keep them for now

        Arguments:
            raw_results (dict): Raw results parsed from self.read_raw_results
            index (str or int): Index or slice of images

        Returns:
            List[results], List[Atoms]

        """
        aimd_results = raw_results.get("aimd", [])
        calc_results = []
        if len(aimd_results) == 0:
            warn(
                "Aimd file is empty! "
                "This is not an error if the calculation "
                "is continued from restart. "
            )
            return None, None

        if isinstance(index, int):
            _images = [aimd_results[index]]
        elif isinstance(index, str):
            _images = aimd_results[string2index(index)]

        ase_images = []
        for result in _images:
            partial_result = {}
            atoms = self.init_atoms.copy()
            if "total energy per atom" in result:
                partial_result["energy"] = result["total energy per atom"] * len(atoms)
            if "free energy per atom" in result:
                partial_result["free energy"] = result["free energy per atom"] * len(
                    atoms
                )

            if "forces" in result:
                # The forces are already re-sorted!
                partial_result["forces"] = result["forces"][self.resort]

            # Modify the atoms in-place
            if "positions" not in result:
                raise ValueError("Cannot have aimd without positions information!")

            atoms.set_positions(
                result["positions"][self.resort], apply_constraint=False
            )

            if "velocities" in result:
                atoms.set_velocities(result["velocities"][self.resort])

            ase_images.append(atoms)
            calc_results.append(partial_result)
        return calc_results, ase_images

    @property
    def sort(self):
        """Wrap the self.sorting dict. If sorting information does not exist,
        use the default slicing
        """

        if self.sorting is None:
            return slice(None, None, None)
        sort = self.sorting.get("sort", [])
        if len(sort) > 0:
            return sort
        else:
            return slice(None, None, None)

    @property
    def resort(self):
        """Wrap the self.sorting dict. If sorting information does not exist,
        use the default slicing
        """

        if self.sorting is None:
            return slice(None, None, None)
        resort = self.sorting.get("resort", [])
        if len(resort) > 0:
            return resort
        else:
            return slice(None, None, None)

    def read_psp_info(self):
        """Parse the psp information from inpt file options
        The psp file locations are relative to the bundle.

        If the files cannot be found, the dict will only contain
        the path
        """
        inpt = self.init_inputs.get("ion", {})
        blocks = inpt.get("atom_blocks", [])
        psp_info = {}
        for block in blocks:
            element = block["ATOM_TYPE"]
            pseudo_path = block["PSEUDO_POT"]
            real_path = (self.directory / pseudo_path).resolve()
            psp_info[element] = {"rel_path": pseudo_path}
            if not real_path.is_file():
                warn(f"Cannot locate pseudopotential {pseudo_path}. ")
            else:
                header = open(real_path, "r").read()
                psp_data = parse_psp8_header(header)
                psp_info[element].update(psp_data)
        return psp_info


def read_sparc(filename, index=-1, include_all_files=True, **kwargs):
    """Parse a SPARC bundle, return an Atoms object or list of Atoms (image)
    with embedded calculator result.

    Arguments:
        filename (str or PosixPath): Filename to the sparc bundle
        index (int or str): Index or slice of the images, following the ase.io.read convention
        include_all_files (bool): If true, parse all output files with indexed suffices
        **kwargs: Additional parameters

    Returns:
       Atoms or List[Atoms]

    """
    # We rely on minimal api version choose, i.e. default or set from env
    api = locate_api()
    sb = SparcBundle(directory=filename, validator=api)
    atoms_or_images = sb.convert_to_ase(
        index=index, include_all_files=include_all_files, **kwargs
    )
    return atoms_or_images


def write_sparc(filename, images, **kwargs):
    """Write sparc file. Images can only be Atoms object
    or list of length 1

    Arguments:
        filename (str or PosixPath): Filename to the output sparc directory
        images (Atoms or List(Atoms)): Atoms object to be written. Only supports writting 1 Atoms
        **kwargs: Additional parameters
    """
    if isinstance(images, Atoms):
        atoms = images
    elif isinstance(images, list):
        if len(images) > 1:
            raise ValueError("SPARC format only supports writing one atoms object!")
        atoms = images[0]
    api = locate_api()
    sb = SparcBundle(directory=filename, mode="w", validator=api)
    sb._write_ion_and_inpt(atoms, **kwargs)
    return


@deprecated(
    "Reading individual .ion file is not recommended. Please use read_sparc instead."
)
def read_sparc_ion(filename, **kwargs):
    """Parse an .ion file inside the SPARC bundle using a wrapper around SparcBundle
    The reader works only when other files (.inpt) exist.

    The returned Atoms object of read_ion method only contains the initial positions

    Arguments:
        filename (str or PosixPath): Filename to the .ion file
        index (int or str): Index or slice of the images, following the ase.io.read convention
        **kwargs: Additional parameters

    Returns:
       Atoms or List[Atoms]
    """
    api = locate_api()
    parent_dir = Path(filename).parent
    sb = SparcBundle(directory=parent_dir, validator=api)
    atoms = sb._read_ion_and_inpt()
    return atoms


# Backward compatibity
read_ion = read_sparc_ion


@deprecated(
    "Writing individual .ion file is not recommended. Please use write_sparc instead."
)
def write_sparc_ion(filename, atoms, **kwargs):
    """Write .ion file using the SparcBundle wrapper. This method will also create the .inpt file

    This is only for backward compatibility

    Arguments:
        filename (str or PosixPath): Filename to the .ion file
        atoms (Atoms): atoms to be written
        **kwargs: Additional parameters
    """
    label = Path(filename).with_suffix("").name
    parent_dir = Path(filename).parent
    api = locate_api()
    sb = SparcBundle(directory=parent_dir, label=label, mode="w", validator=api)
    sb._write_ion_and_inpt(atoms, **kwargs)
    return atoms


# Backward compatibility
write_ion = write_sparc_ion


@deprecated(
    "Reading individual .static file is not recommended. Please use read_sparc instead."
)
def read_sparc_static(filename, index=-1, **kwargs):
    """Parse a .static file bundle using a wrapper around SparcBundle
    The reader works only when other files (.ion, .inpt) exist.

    Arguments:
        filename (str or PosixPath): Filename to the .static file
        index (int or str): Index or slice of the images, following the ase.io.read convention
        **kwargs: Additional parameters

    Returns:
       Atoms or List[Atoms]
    """
    parent_dir = Path(filename).parent
    api = locate_api()
    sb = SparcBundle(directory=parent_dir, validator=api)
    # In most of the cases the user wants to inspect all images
    kwargs = kwargs.copy()
    if "include_all_files" not in kwargs:
        kwargs.update(include_all_files=True)
    atoms_or_images = sb.convert_to_ase(index=index, **kwargs)
    return atoms_or_images


# Backward compatibility
read_static = read_sparc_static


@deprecated(
    "Reading individual .geopt file is not recommended. Please use read_sparc instead."
)
def read_sparc_geopt(filename, index=-1, **kwargs):
    """Parse a .geopt file bundle using a wrapper around SparcBundle
    The reader works only when other files (.ion, .inpt) exist.

    Arguments:
        filename (str or PosixPath): Filename to the .geopt file
        index (int or str): Index or slice of the images, following the ase.io.read convention
        **kwargs: Additional parameters

    Returns:
       Atoms or List[Atoms]
    """
    parent_dir = Path(filename).parent
    api = locate_api()
    sb = SparcBundle(directory=parent_dir, validator=api)
    kwargs = kwargs.copy()
    if "include_all_files" not in kwargs:
        kwargs.update(include_all_files=True)
    atoms_or_images = sb.convert_to_ase(index=index, **kwargs)
    return atoms_or_images


# Backward compatibility
read_geopt = read_sparc_geopt


@deprecated(
    "Reading individual .aimd file is not recommended. Please use read_sparc instead."
)
def read_sparc_aimd(filename, index=-1, **kwargs):
    """Parse a .static file bundle using a wrapper around SparcBundle
    The reader works only when other files (.ion, .inpt) exist.

    Arguments:
        filename (str or PosixPath): Filename to the .aimd file
        index (int or str): Index or slice of the images, following the ase.io.read convention
        **kwargs: Additional parameters

    Returns:
       Atoms or List[Atoms]
    """
    parent_dir = Path(filename).parent
    api = locate_api()
    sb = SparcBundle(directory=parent_dir, validator=api)
    kwargs = kwargs.copy()
    if "include_all_files" not in kwargs:
        kwargs.update(include_all_files=True)
    atoms_or_images = sb.convert_to_ase(index=index, **kwargs)
    return atoms_or_images


# Backward compatibility
read_aimd = read_sparc_aimd


def __register_new_filetype():
    """Register the filetype() function that allows recognizing .sparc as directory
    This method should only be called for ase==3.22 compatibility and for ase-gui
    In future versions of ase gui where format is supported, this method should be removed
    """
    import sys

    from ase.io import formats as hacked_formats
    from ase.io.formats import filetype as _old_filetype
    from ase.io.formats import ioformats

    def _new_filetype(filename, read=True, guess=True):
        """A hacked solution for the auto format recovery"""
        path = Path(filename)
        ext = path.name
        if ".sparc" in ext:
            return "sparc"
        else:
            if path.is_dir():
                if (len(list(path.glob("*.ion"))) > 0) and (
                    len(list(path.glob("*.inpt"))) > 0
                ):
                    return "sparc"
            return _old_filetype(filename, read, guess)

    hacked_formats.filetype = _new_filetype
    sys.modules["ase.io.formats"] = hacked_formats
    return


@deprecated(
    "register_ase_io_sparc will be deprecated for future releases. Please upgrade ase>=3.23."
)
def register_ase_io_sparc(name="sparc"):
    """
    **Legacy register of io-formats for ase==3.22**
    **For ase>=3.23, use the package entrypoint registration**
    Monkey patching the ase.io and ase.io.formats
    So that the following formats can be used
    after `import sparc`

    ```
    from ase.io import sparc
    ase.io.read("test.sparc")
    atoms.write("test.sparc")
    ```

    The register method only aims to work for ase 3.22
    the develope version of ase provides a much more powerful
    register mechanism, we can wait.
    """
    import sys
    from warnings import warn

    import pkg_resources
    from ase.io.formats import define_io_format as F
    from ase.io.formats import ioformats

    name = name.lower()
    if name in ioformats.keys():
        return
    desc = "SPARC .sparc bundle"

    # Step 1: patch the ase.io.sparc module
    try:
        entry_points = next(
            ep for ep in pkg_resources.iter_entry_points("ase.io") if ep.name == "sparc"
        )
        _monkey_mod = entry_points.load()
    except Exception as e:
        warn(
            (
                "Failed to load entrypoint `ase.io.sparc`, "
                "you may need to reinstall sparc python api.\n"
                "You may still use `sparc.read_sparc` and "
                "`sparc.write_sparc` methods, "
                "but not `ase.io.read`\n",
                f"The error is {e}",
            )
        )
        return

    sys.modules[f"ase.io.{name}"] = _monkey_mod
    __register_new_filetype()

    # Step 2: define a new format
    F(
        name,
        desc=desc,
        code="+S",  # read_sparc has multi-image support
        ext="sparc",
    )

    if name not in ioformats.keys():
        warn(
            (
                "Registering .sparc format with ase.io failed. "
                "You may still use `sparc.read_sparc` and "
                "`sparc.write_sparc` methods. \n"
                "Please contact the developer to report this issue."
            )
        )
        return

    import tempfile

    from ase.io import read

    with tempfile.TemporaryDirectory(suffix=".sparc") as tmpdir:
        try:
            read(tmpdir.name)
        except Exception as e:
            emsg = str(e).lower()
            if "bundletrajectory" in emsg:
                warn(
                    "Atomatic format inference for sparc is not correctly registered. "
                    "You may need to use format=sparc in ase.io.read and ase.io.write. "
                )
    # Add additional formats including .ion (r/w), .static, .geopt, .aimd
    F(
        "ion",
        desc="SPARC .ion file",
        module="sparc",
        code="1S",
        ext="ion",
    )
    F(
        "static",
        desc="SPARC single point results",
        module="sparc",
        code="+S",
        ext="static",
    )
    F(
        "geopt",
        desc="SPARC geometric optimization results",
        module="sparc",
        code="+S",
        ext="geopt",
    )
    F("aimd", desc="SPARC AIMD results", module="sparc", code="+S", ext="aimd")

    # TODO: remove print options as it may be redundant
    print("Successfully registered sparc formats with ase.io!")


# ase>=3.23 uses new ExternalIOFormat as registered entrypoints
# Please do not use from ase.io.formats import ExternalIOFormat!
# This causes circular import
try:
    from ase.utils.plugins import ExternalIOFormat as EIF
except ImportError:
    # Backward Compatibility
    from typing import List, NamedTuple, Optional, Union

    # Copy definition from 3.23
    # Name is defined in the entry point
    class ExternalIOFormat(NamedTuple):
        desc: str
        code: str
        module: Optional[str] = None
        glob: Optional[Union[str, List[str]]] = None
        ext: Optional[Union[str, List[str]]] = None
        magic: Optional[Union[bytes, List[bytes]]] = None
        magic_regex: Optional[bytes] = None

    EIF = ExternalIOFormat

format_sparc = EIF(
    desc="SPARC .sparc bundle",
    module="sparc.io",
    code="+S",  # read_sparc has multi-image support
    ext="sparc",
)
format_ion = EIF(
    desc="SPARC .ion file",
    module="sparc.io",
    code="1S",
    ext="ion",
)
format_static = EIF(
    desc="SPARC single point results",
    module="sparc.io",
    code="+S",
    glob=["*.static", "*.static_*"],
)
format_geopt = EIF(
    desc="SPARC geometric optimization results",
    module="sparc.io",
    code="+S",
    glob=["*.geopt", "*.geopt_*"],
)
format_aimd = EIF(
    desc="SPARC AIMD results",
    module="sparc",
    code="+S",
    glob=["*.aimd*", "*.geopt_*"],
)
