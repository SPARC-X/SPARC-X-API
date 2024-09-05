# Temporary fix to recognize .sparc from commandline
from .io import __register_new_filetype

__register_new_filetype()


# The cli part should be considered deprecated now.
def main(
    prog="sparc-ase",
    description="ASE command line tool with SPARC support",
    hook=None,
    args=None,
):
    import sys

    from ase.cli.main import main as ase_main

    ase_main(prog=prog, description=description, hook=hook, args=args)


if __name__ == "__main__":
    main()
