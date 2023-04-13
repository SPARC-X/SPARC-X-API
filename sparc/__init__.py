# __Author__ = "Ben Comer"
# from .sparc_core import SPARC
from .sparc_io_bundle import read_sparc, write_sparc
from .sparc_io_bundle import register_ase_io_sparc

register_ase_io_sparc()
