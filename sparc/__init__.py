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
    from packaging import version

    from .calculator import SPARC
    from .io import read_sparc, register_ase_io_sparc, write_sparc

    # If ase version less than 3.23, use manual register function
    # Otherwise use the new entry point
    if version.parse(ase.__version__) < version.parse("3.23"):
        register_ase_io_sparc()
    else:
        # register calculator class <experimental>
        from ase.calculators.calculator import register_calculator_class

        register_calculator_class("sparc", SPARC)
else:
    # If importing is not complete, any code trying to directly import
    # the following attributes will raise ImportError
    read_sparc = _missing_deps_func
    write_sparc = _missing_deps_func
    SPARC = SPARCMissingDeps
