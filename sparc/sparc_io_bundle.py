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
from .sparc_parsers.pseudopotential import copy_psp_file
from .common import psp_dir as default_psp_dir
from .download_data import is_psp_download_complete


class SparcBundle:
    """Provide access to a calculation folder of SPARC as a simple bundle

    The bundle can be optionally named as .sparc following the ASE's .bundle format

    Currently the write method only supports 1 image, while read method support reading
    atoms results in following conditions

    1) No calculation (minimal): .ion + .inpt file --> 1 image
    2) Single point calculation: .ion + .inpt + .out + .static --> 1 image with calc
    3) Multiple SP calculations: chain all .out{digits} and .static{digitis} outputs
    4) Relaxation: read from .geopt and .out (supporting chaining)
    5) AIMD: read from .aimd and .out (support chaining)


    Currently, the bundle object is intended to be used for one-time read / write

    TODO: multiple occurance support
    TODO: archive support



    """

    psp_env = ["SPARC_PSP_PATH", "SPARC_PP_PATH"]

    def __init__(self, directory, mode="r", atoms=None, label=None, psp_dir=None):
        self.directory = Path(directory)
        # TODO: more sensible naming for name?
        self.prefix = self.directory.resolve().with_suffix("").name
        self.label = self._make_label(label)  # name of the main sparc file
        self.mode = mode.lower()
        assert self.mode in (
            "r",
            "w",
            "a",
        ), f"Invalid mode {self.mode}! Must one of 'r', 'w' or 'a'"
        # TODO: assigning atoms here is probably not useful!
        self.atoms = atoms
        self.psp_dir = self.__find_psp_dir(psp_dir)
        # Sorting should be consistent across the whole bundle!
        self.sorting = None
        self.last_image = -1

    def _find_files(self):
        """Find all files matching '{label}.*'"""
        return list(self.directory.glob(f"{self.label}.*"))

    def _make_label(self, label=None):
        illegal_chars = '\\/:*?"<>|'
        label_ = label if label is not None else self.prefix
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
                        "PSP directory bundled with sparc-dft-api is broken! Please use `sparc.download_data` to re-download them!"
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
        """Read the ion and inpt files together"""
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
        **kwargs,
    ):
        """Write the ion and inpt files to a bundle. This method only supports writing 1 image.
        If input_parameters are empty, there will only be .ion writing the positions and .inpt writing a minimal cell information

        """
        if self.mode != "w":
            raise ValueError(
                "Cannot write input files while sparc bundle is opened in read or append mode!"
            )
        os.makedirs(self.directory, exist_ok=True)
        atoms = self.atoms.copy() if atoms is None else atoms.copy()
        # TODO: make the parameter more explicit
        pseudopotentials = kwargs.pop("pseudopotentials", {})
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

    def read_results(self, image=-1):
        """Parse all files using the given self.label.
        The results are merged dict from all file formats

        For now, we just read the FIRST appearing files to confirm it's working

        TODO: support multi occurance

        """
        results_dict = {}
        # Find the max output index
        # TODO: move this into another function
        last_out = sorted(self.directory.glob(f"{self.label}.out*"), reverse=True)[0]
        print("Last output file: ", last_out)
        suffix = last_out.suffix
        if suffix == ".out":
            self.last_image = 0
        else:
            self.last_image = int(suffix.split("_")[1])
        current_image = range(self.last_image + 1)[image]

        print("Current image to read: ", current_image)
        for ext in ("ion", "inpt"):
            # TODO: maybe a hack for adding multiple occruance here
            f = self._indir(ext, occur=0)
            if f.is_file():
                data_dict = globals()[f"_read_{ext}"](f)
                results_dict.update(data_dict)
        for ext in ("geopt", "static", "aimd", "out"):
            f = self._indir(ext, occur=current_image)
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
                assert (tuple(self.sorting["sort"]) == tuple(sorting["sort"])) and (
                    tuple(self.sorting["resort"]) == tuple(sorting["resort"])
                ), "Sorting information changed!"
        return results_dict


def read_sparc(filename, *args, **kwargs):
    """Very simple PoC now"""
    sb = SparcBundle(directory=filename)
    atoms = sb._read_ion_and_inpt()
    return atoms


def write_sparc(filename, atoms, **kwargs):
    sb = SparcBundle(directory=filename)
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

    The register method only aims to work for ase <= 3.22
    the develope version of ase provides a much more powerful
    register mechanism, we can wait.
    """
    from ase.io.formats import define_io_format as F
    from ase.io.formats import ioformats
    from importlib import import_module
    import pkg_resources
    import sys
    from warnings import warn

    name = name.lower()
    if name in ioformats.keys():
        return
    desc = "Bundled calculation directory for SPARC " "quantum chemistry code"

    # Step 1: patch the ase.io.sparc module
    try:
        entry_points = next(
            ep for ep in pkg_resources.iter_entry_points("ase.io") if ep.name == "sparc"
        )
        _monkey_mod = entry_points.load()
    except Exception as e:
        warn(
            (
                "Failed to load entrypoint `ase.io.sparc`, you may need to reinstall sparc python api.\n"
                'You may still use `sparc.read_sparc` and `sparc.write_sparc` methods, but not `ase.io.read("test.sparc")`\n',
                f"The error is {e}",
            )
        )
        return

    sys.modules[f"ase.io.{name}"] = _monkey_mod

    # Step 2: define a new format
    F(
        name,
        desc=desc,
        code="1S",  # Currently make only 1 image
        ext="sparc",
    )

    if name not in ioformats.keys():
        warn(
            (
                "Registering .sparc format with ase.io failed. "
                "You may still use `sparc.read_sparc` and `sparc.write_sparc` methods. "
                "You're welcome to contact the developer to report this issue."
            )
        )
        return

    # TODO: remove print options as it may be redundant
    print("Successfully registered sparc format with ase.io!")
