import pkg_resources
from pathlib import Path



print(data_dir, repo_dir)


from .sparc_io_bundle import read_sparc, write_sparc
from .sparc_io_bundle import register_ase_io_sparc

register_ase_io_sparc()
