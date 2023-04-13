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


# @reader
def read_sparc(filename, *args, **kwargs):
    # raise NotImplementedError
    print("I'm the real writer")
    pass


# @writer
def write_sparc(filename, atoms, label=None, sort=None, copy_psp=True, **kwargs):
    # raise NotImplementedError
    print("I'm the real writer")
    pass


####################################################
# Monkey patching the ase.io and ase.io.formats
# So that the following formats can be used
# after `import sparc`
# ase.io.read("test.sparc")
# atoms.write("test.sparc")
# from ase.io import sparc
####################################################
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
