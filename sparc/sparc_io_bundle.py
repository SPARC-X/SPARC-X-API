"""Providing a new bundled SPARC file format
.sparc

Many of the logics are taken from
ase.io.vasp

ase.io.trajectory

"""
import re

import numpy as np

from ase import Atoms, Atom
from ase.utils import reader, writer
from ase.io.utils import ImageIterator
from ase.io import ParseError
from pathlib import Path

# from .sparc_parsers.ion import read_ion, write_ion

from ase import io
import numpy as np
from pathlib import Path
import shutil
import tarfile

from .sparc_parsers.ion import _read_ion, _write_ion, _ion_coord_to_ase_pos
from .sparc_parsers.inpt import _read_inpt, _write_inpt, _inpt_cell_to_ase_cell
from .sparc_parsers.atoms import dict_to_atoms, atoms_to_dict
from .inputs import SparcInputs

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

    """

    def __init__(self, directory, mode="r", atoms=None, label=None):
        self.directory = Path(directory)
        # TODO: more sensible naming for name?
        self.name = self.directory.with_suffix("").name
        self.label = label if label is not None else self.name
        self.file_list = []  # list of file names in the directory
        self.sparc_file = None  # name of the main sparc file
        self.mode = mode.lower()
        self.atoms = atoms

    def _indir(self, ext, label=None):
        """Find the file with {label}.{ext} under current dir

        if label is None, use the default
        # TODO: how about recursive?
        """
        label = self.label if label is None else label
        if not ext.startswith("."):
            ext = "." + ext
        target = self.directory / f"{label}{ext}"
        return target

    def _read_ion_and_inpt(self):
        """Read the ion and inpt files together"""
        f_ion, f_inpt = self._indir(".ion"), self._indir(".inpt")
        ion_data = _read_ion(f_ion)
        inpt_data = _read_inpt(f_inpt)
        merged_data = {**ion_data, **inpt_data}
        return dict_to_atoms(merged_data)

    def _write_ion_and_inpt(self, atoms=None, label=None, direct=False, sort=None, ignore_constraints=False, wrap=False,
                            # Below are the parameters from v1
                            # scaled -> direct, ignore_constraints --> not add_constraints
                             scaled=False, add_constraints=True, copy_psp=True, comment="", input_parameters={}):
        """Write the ion and inpt files to a bundle. This method only supports writing 1 image.
        If input_parameters are empty, there will only be .ion writing the positions and .inpt writing a minimal cell information

        """
        atoms = self.atoms.copy() if atoms is None else atoms.copy()
        
        

    # def _determine_state(self):
    #     """Determine the state of files inside the bundle according to the existence of files

    #     1. .ion / .inpt only --> initial structure (no calc)
    #     2. existing output file --> reading output has higher priority 
    #     """


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
