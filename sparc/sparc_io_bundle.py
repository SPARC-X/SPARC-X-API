"""Providing a new bundled SPARC file format
.sparc

Many of the logics are taken from
ase.io.vasp

ase.io.trajectory

"""
import re

import numpy as np

from ase import Atoms
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
        self.directory_path = Path(directory)
        # TODO: more sensible naming for name?
        self.name = self.directory_path.with_suffix("").name
        self.label = label if label is not None else self.name
        self.file_list = []  # list of file names in the directory
        self.sparc_file = None  # name of the main sparc file
        self.mode = mode.lower()
        self.atoms = atoms

        # if mode == 'r':
        #     if self.directory_path.is_file() and self.directory_path.suffix == '.tar':
        #         self._extract_tar()
        #         self.directory_path = self.directory_path.with_suffix('')
        #     self._validate_directory()
        # elif mode == 'w':
        #     self.directory_path.mkdir(parents=True, exist_ok=True)

        # if atoms is not None:
        #     self.write_atoms(atoms)

    def read_sparc(self):
        # read the main sparc file and return the data
        pass

    def write_sparc(self, data, copy_psp=None):
        # write the data to the main sparc file
        pass

    def read_atoms(self, filename):
        # read atoms from a file
        return io.read(filename)

    def write_atoms(self, atoms, filename=None, copy_psp=None):
        # write atoms to a file
        if filename is None:
            if self.label is None:
                filename = self.directory_path / "atoms.xyz"
            else:
                filename = self.directory_path / f"{self.label}.xyz"
        io.write(str(filename), atoms)
        self.file_list.append(filename.name)

        if copy_psp is not None:
            for psp_path in copy_psp:
                if not isinstance(psp_path, Path):
                    psp_path = Path(psp_path)

                if not psp_path.exists():
                    raise ValueError(f"{psp_path} does not exist")

                psp_filename = psp_path.name
                shutil.copy(str(psp_path), str(self.directory_path))
                if not (self.directory_path / psp_filename).exists():
                    raise ValueError(
                        f"Failed to copy {psp_filename} to {self.directory_path}"
                    )
                self.file_list.append(psp_filename)

    def _validate_directory(self):
        if not self.directory_path.exists():
            raise ValueError(f"Directory {self.directory_path} does not exist")

        sparc_files = []
        for file_name in self.directory_path.iterdir():
            if file_name.suffix == ".sparc":
                sparc_files.append(file_name)
        if not sparc_files:
            raise ValueError("No sparc files found in the directory")
        elif len(sparc_files) > 1:
            raise ValueError("More than one sparc files found in the directory")

        self.sparc_file = sparc_files[0]

        if self.label is not None:
            label_files = []
            for file_name in self.directory_path.iterdir():
                if (
                    file_name.name.startswith(f"{self.label}.")
                    and not file_name.suffix == ".sparc"
                ):
                    label_files.append(file_name)
            if not label_files:
                raise ValueError(f"No {self.label} files found in the directory")

        self.file_list = [file.name for file in self.directory_path.iterdir()]

    def validate(self):
        pass
        # validate if the directory contains minimal files to be a valid sparc


def read_sparc(filename, *args, **kwargs):
    pass


def write_sparc(filename, atoms, label=None, sort=None, copy_psp=True, **kwargs):
    pass


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
