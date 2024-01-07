"""Initialization of sparc-x-api

For submodules like download_data and api, ase / numpy may be ignored,
and run using standard python libaries. This may be useful for cases like
conda build and CI where not all dependencies are present
"""


def _missing_deps_func(*args, **kwargs):
    raise ImportError("Importing fails for ase / numpy!")


class SPARCMissingDeps:
    def __init__(self, *args, **kwargs):
        raise ImportError(
            "Cannot initialize sparc.SPARC because the required dependencies (ase and numpy) are not available."
        )

    def __getattr__(self, name):
        raise ImportError(
            f"Cannot access '{name}' on sparc.SPARC because the required dependencies (ase and numpy) are not available."
        )


try:
    import ase
    import numpy

    _import_complete = True
except ImportError:
    _import_complete = False

if _import_complete:
    from .calculator import SPARC
    from .io import read_sparc, register_ase_io_sparc, write_sparc

    register_ase_io_sparc()
else:
    # If importing is not complete, any code trying to directly import
    # the following attributes will raise ImportError
    read_sparc = _missing_deps_func
    write_sparc = _missing_deps_func
    SPARC = SPARCMissingDeps
