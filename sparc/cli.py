def main(
    prog="sparc-ase",
    description="ASE command line tool with SPARC support",
    hook=None,
    args=None,
):
    import sys

    from ase.cli.main import main as ase_main
    from ase.io import sparc as mod_sparc
    from ase.io.sparc import read_sparc as _old_read_sparc

    import sparc

    def _new_read_sparc(filename, index, **kwargs):
        return _old_read_sparc(filename, index, include_all_files=True, **kwargs)

    try:
        sys.modules["ase.io.sparc"].read_sparc = _new_read_sparc
        ase_main(prog=prog, description=description, hook=hook, args=args)
    finally:
        sys.modules["ase.io.sparc"].read_sparc = _old_read_sparc


if __name__ == "__main__":
    main()
