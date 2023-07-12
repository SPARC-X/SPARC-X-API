"""Providing a new bundled SPARC file format
.sparc

Many of the logics are taken from
ase.io.vasp

ase.io.trajectory

"""
import os

from pathlib import Path

# from .sparc_parsers.ion import read_ion, write_ion

from warnings import warn

# various io formatters
from .sparc_parsers.ion import _read_ion, _write_ion
from .sparc_parsers.inpt import _read_inpt, _write_inpt
from .sparc_parsers.static import _read_static
from .sparc_parsers.geopt import _read_geopt
from .sparc_parsers.aimd import _read_aimd
from .sparc_parsers.out import _read_out

from .sparc_parsers.atoms import dict_to_atoms, atoms_to_dict
from .sparc_parsers.pseudopotential import copy_psp_file, parse_psp8_header
from .common import psp_dir as default_psp_dir
from .download_data import is_psp_download_complete
from .utils import string2index

from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.units import Hartree
from ase.atoms import Atoms


class SparcBundle:
    """Provide access to a calculation folder of SPARC as a simple bundle

    The bundle can be optionally named as .sparc following the ASE's
    .bundle format

    Currently the write method only supports 1 image, while read method support reading
    atoms results in following conditions

    1) No calculation (minimal): .ion + .inpt file --> 1 image 2)
    Single point calculation: .ion + .inpt + .out + .static --> 1
    image with calc 3) Multiple SP calculations: chain all
    .out{digits} and .static{digitis} outputs 4) Relaxation: read from
    .geopt and .out (supporting chaining) 5) AIMD: read from .aimd and
    .out (support chaining)


    Currently, the bundle object is intended to be used for one-time
    read / write

    TODO: multiple occurance support
    TODO: archive support

    """

    psp_env = ["SPARC_PSP_PATH", "SPARC_PP_PATH"]

    def __init__(
        self, directory, mode="r", atoms=None, label=None, psp_dir=None,
    ):
        self.directory = Path(directory)
        self.mode = mode.lower()
        assert self.mode in (
            "r",
            "w",
            "a",
        ), f"Invalid mode {self.mode}! Must one of 'r', 'w' or 'a'"
        self.label = self._make_label(label)  # name of the main sparc file
        # TODO: assigning atoms here is probably not useful!
        self.init_atoms = atoms.copy() if atoms is not None else None
        self.init_inputs = {}
        self.psp_data = {}
        self.raw_results = {}
        self.psp_dir = self.__find_psp_dir(psp_dir)
        # Sorting should be consistent across the whole bundle!
        self.sorting = None
        self.last_image = -1
        self.allow_custom_parameters = allow_custom_parameters

    def _find_files(self):
        """Find all files matching '{label}.*'"""
        return list(self.directory.glob(f"{self.label}.*"))

    def _make_label(self, label=None):
        """Infer the label from the bundle

        Special cases if label is None:
        1. read mode --> get the ion file name
        2. write mode --> infer from the directory
        """
        # TODO: more sensible naming for name?
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
                # TODO: customize error msg
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
                        "PSP directory bundled with sparc-dft-api is broken! "
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
        """Find the file with {label}.{ext} under current dir

        if label is None, use the default
        # TODO: how about recursive?
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
        """Read the ion and inpt files together

        This method should be rarely used
        """
        f_ion, f_inpt = self._indir(".ion"), self._indir(".inpt")
        ion_data = _read_ion(f_ion)
        inpt_data = _read_inpt(f_inpt)
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
        custom_parameters={},
        **kwargs,
    ):
        """Write the ion and inpt files to a bundle. This method only
        supports writing 1 image.  If input_parameters are empty,
        there will only be .ion writing the positions and .inpt
        writing a minimal cell information

        """
        if self.mode != "w":
            raise ValueError(
                "Cannot write input files while sparc bundle is opened in read or append mode!"
            )
        os.makedirs(self.directory, exist_ok=True)
        atoms = self.atoms.copy() if atoms is None else atoms.copy()
        # TODO: make the parameter more explicit
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
        # TODO: special input param handling
        data_dict["inpt"]["params"].update(merged_inputs)
        # TODO: label

        # If copy_psp, change the PSEUDO_POT field and copy the files
        if copy_psp:
            for block in data_dict["ion"]["atom_blocks"]:
                if "PSEUDO_POT" in block:
                    origin_psp = block["PSEUDO_POT"]
                    target_dir = self.directory
                    target_fname = copy_psp_file(origin_psp, target_dir)
                    block["PSEUDO_POT"] = target_fname

        _write_ion(self._indir(".ion"), data_dict)
        _write_inpt(self._indir(".inpt"), data_dict)
        return

    def read_raw_results(self, include_all_files=False):
        """Parse all files using the given self.label.
        The results are merged dict from all file formats

        Argument all_files: True --> include all files (out, out_01,
        out_02, etc) when all files are included, output is a list of
        dicts; otherwise a single dict

        """
        # Find the max output index
        # TODO: move this into another function
        last_out = sorted(
            self.directory.glob(f"{self.label}.out*"), reverse=True
        )
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

        # print(self.last_image, self.num_calculations)

        if include_all_files:
            results = [
                self._read_results_from_index(index)
                for index in range(self.num_calculations)
            ]
        else:
            results = self._read_results_from_index(self.last_image)

        self.raw_results = results

        if include_all_files:
            init_raw_results = self.raw_results[0]
            # self.sorting = self.raw_results[0]["ion"]["sorting"]
        else:
            init_raw_results = self.raw_results.copy()
            # self.sorting = self.raw_results["ion"]["sorting"]

        # TODO: init is actually last!
        self.init_atoms = dict_to_atoms(init_raw_results)
        self.init_inputs = {
            "ion": init_raw_results["ion"],
            "inpt": init_raw_results["inpt"],
        }
        self.psp_data = self.read_psp_info()
        return self.raw_results

    def _read_results_from_index(self, index, d_format="{:02d}"):
        """Read the results from one calculation index, and return a
        single raw result dict"""
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
        # TODO: make another function to check sanity
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
                assert (
                    tuple(self.sorting["sort"]) == tuple(sorting["sort"])
                ) and (
                    tuple(self.sorting["resort"]) == tuple(sorting["resort"])
                ), "Sorting information changed!"
        return results_dict

    def convert_to_ase(self, index=-1, include_all_files=False, **kwargs):
        """Read the raw results from the bundle and create atoms with
        single point calculators

        TODO: what to do about the indices?

        """
        # Convert to images!
        rs = self.read_raw_results(include_all_files=include_all_files)
        if isinstance(rs, dict):
            raw_results = [rs]
        else:
            raw_results = list(rs)
        res_images = []
        # print("RAW RES: ", raw_results)
        for entry in raw_results:
            if "static" in entry:
                calc_results, images = self._extract_static_results(
                    entry, index=":"
                )
            elif "geopt" in entry:
                calc_results, images = self._extract_geopt_results(
                    entry, index=":"
                )
            elif "aimd" in entry:
                calc_results, images = self._extract_aimd_results(
                    entry, index=":"
                )
            else:
                calc_results, images = None, [self.init_atoms.copy()]

            if images is not None:
                if calc_results is not None:
                    images = self._make_singlepoint(
                        calc_results, images, entry)
                res_images.extend(images)

        if isinstance(index, int):
            return res_images[index]
        else:
            return res_images[string2index(index)]

    def _make_singlepoint(self, calc_results, images, raw_results):
        """Convert a calculator dict and images of Atoms to list of
        SinglePointDFTCalculators

        The calculator also takes parameters from ion, inpt that exist
        in self.raw_results

        """
        converted_images = []
        for res, _atoms in zip(calc_results, images):
            atoms = _atoms.copy()
            sp = SinglePointDFTCalculator(atoms)
            # Simply copy the results?
            sp.results.update(res)
            sp.name = "sparc"
            sp.kpts = (
                raw_results["inpt"]
                .get("parameters", {})
                .get("KPOINT_GRID", None)
            )
            # There may be a better way handling the parameters...
            sp.parameters = raw_results["inpt"].get("parameters", {})
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

        """
        # TODO: implement the multi-file static
        static_results = raw_results.get("static", {})
        calc_results = {}
        if "free energy" in static_results:
            calc_results["energy"] = static_results["free energy"]
            calc_results["free energy"] = static_results["free energy"]

        if "forces" in static_results:
            # TODO: what about non-sorted ones
            calc_results["forces"] = static_results["forces"][self.resort]

        if "stress" in static_results:
            calc_results["stress"] = static_results["stress"]

        atoms = self.init_atoms.copy()
        if "atoms" in static_results:
            atoms_dict = static_results["atoms"]
            # TODO: detect change in atomic symbols!
            # TODO: Check naming, is it coord_frac or scaled_positions?
            if "coord_frac" in atoms_dict:
                # TODO: check if set_scaled_positions requires constraint?
                atoms.set_scaled_positions(
                    atoms_dict["coord_frac"][self.resort]
                )
            elif "coord" in atoms_dict:
                atoms.set_positions(
                    atoms_dict["coord"][self.resort], apply_constraint=False
                )
        return [calc_results], [atoms]

    def _extract_geopt_results(self, raw_results, index=":"):
        """Extract the static calculation results and atomic
        structure(s) Returns: calc_results: dict with at least energy
        value atoms: ASE atoms object The priority is to parse
        position from static file first, then fallback from ion + inpt

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
        # import pdb; pdb.set_trace()
        for result in _images:
            atoms = self.init_atoms.copy()
            partial_result = {}
            if "energy" in result:
                partial_result["energy"] = result["energy"]
                # TODO: shall we distinguish?
                partial_result["free energy"] = result["energy"]

            if "forces" in result:
                # TODO: what about non-sorted calculations
                partial_result["forces"] = result["forces"][self.resort]

            if "stress" in result:
                partial_result["stress"] = result["stress"]

            # Modify the atoms copy
            if "positions" not in result:
                raise ValueError(
                    "Cannot have geopt without positions information!"
                )
            atoms.set_positions(
                result["positions"][self.resort], apply_constraint=False
            )
            if "ase_cell" in result:
                atoms.set_cell(result["ase_cell"])
            calc_results.append(partial_result)
            ase_images.append(atoms)

        return calc_results, ase_images

    def _extract_aimd_results(self, raw_results, index=":"):
        """Extract energy / forces from aimd results

        For calculator, we only need the last image

        We probably want more information for the AIMD calculations,
        but I'll keep them for now

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
                partial_result["energy"] = result[
                    "total energy per atom"
                ] * len(atoms)
            if "free energy per atom" in result:
                partial_result["free energy"] = result[
                    "free energy per atom"
                ] * len(atoms)

            if "forces" in result:
                # The forces are already re-sorted!
                partial_result["forces"] = result["forces"][self.resort]

            # Modify the atoms in-place
            if "positions" not in result:
                raise ValueError(
                    "Cannot have aimd without positions information!"
                )

            atoms.set_positions(
                result["positions"][self.resort], apply_constraint=False
            )

            # TODO: need to get an example for NPT MD to set Cell
            # TODO: need to set stress information

            if "velocities" in result:
                atoms.set_velocities(result["velocities"][self.resort])

            ase_images.append(atoms)
            calc_results.append(partial_result)
        return calc_results, ase_images

    # def get_ionic_steps(self, raw_results):
    #     """Get last ionic step dict from raw results"""
    #     out_results = raw_results.get("out", {})
    #     ionic_steps = out_results.get("ionic_steps", [])
    #     return ionic_steps

    # def _extract_output_results(self, raw_results):
    #     """Extract extra information from results, need to be more polished
    #     (maybe move to calculator?)
    #     """
    #     last_step = self.get_ionic_step(raw_results)[-1]
    #     if "fermi level" in last_step:
    #         value = last_step["fermi level"]["value"]
    #         unit = last_step["fermi level"]["unit"]
    #         if unit.lower() == "ev":
    #             self.results["fermi"] = value
    #         # Should rarely happen, but keep it here!
    #         elif unit.lower() == "hartree":
    #             self.results["fermi"] = value * Hartree
    #         else:
    #             raise ValueError("Wrong unit in Fermi!")
    #     return

    @property
    def sort(self):
        """wrap the self.sorting dict. If sorting information does not exist,
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
        """wrap the self.sorting dict. If sorting information does not exist,
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


def read_sparc(filename, index=-1, include_all_files=False, **kwargs):
    """Parse a SPARC bundle, return an Atoms object or list of Atoms (image)
    with embedded calculator result.

    """
    sb = SparcBundle(directory=filename)
    atoms_or_images = sb.convert_to_ase(
        index=index, include_all_files=include_all_files, **kwargs
    )
    return atoms_or_images


def write_sparc(filename, images, **kwargs):
    """Write sparc file. Images can only be Atoms object
    or list of length 1
    """
    if isinstance(images, Atoms):
        atoms = images
    elif isinstance(images, list):
        if len(images) > 1:
            raise ValueError(
                "SPARC format only supports writing one atoms object!"
            )
        atoms = images[0]
    sb = SparcBundle(directory=filename, mode="w")
    sb._write_ion_and_inpt(atoms, **kwargs)
    return


def register_ase_io_sparc(name="sparc"):
    """
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
    from ase.io.formats import define_io_format as F
    from ase.io.formats import ioformats
    from ase.io.formats import filetype as _old_filetype
    from ase.io import formats as hacked_formats

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

    import pkg_resources
    import sys
    from warnings import warn

    name = name.lower()
    if name in ioformats.keys():
        return
    desc = "SPARC .sparc bundle"

    # Step 1: patch the ase.io.sparc module
    try:
        entry_points = next(
            ep
            for ep in pkg_resources.iter_entry_points("ase.io")
            if ep.name == "sparc"
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

    hacked_formats.filetype = _new_filetype

    sys.modules[f"ase.io.{name}"] = _monkey_mod
    sys.modules["ase.io.formats"] = hacked_formats
    # sys.modules[f"ase.io.format"] = _monkey_mod

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

    from ase.io import read
    import tempfile

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

    # TODO: remove print options as it may be redundant
    print("Successfully registered sparc format with ase.io!")
