"""A simple test module for sparc python api
Usage:
python -m sparc.quicktest
A few tests will be performed
1) Is the sparc format recognizable by ASE
2) Is the command properly set
3) Is the psp directory accessible
4) Is the json api accessible
5) Can the command actually run a very simple calculation
"""
from .utils import cprint


def import_test():
    cprint("Testing import...", color="COMMENT")
    import sparc
    from ase.io.formats import ioformats

    if "sparc" not in ioformats.keys():
        return False
    try:
        from ase.io import sparc
    except ImportError:
        return False
    return hasattr(sparc, "read_sparc") and hasattr(sparc, "write_sparc")


def psp_test():
    cprint("Testing pseudo potential path...", color="COMMENT")
    from sparc.io import SparcBundle
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        sb = SparcBundle(directory=tmpdir)
        psp_dir = sb.psp_dir
        if psp_dir is not None:
            cprint(f"Found psp files at {psp_dir}", color="OKBLUE")
            return True
        else:
            cprint(
                (
                    "No psp files found! \n"
                    "Please make sure you have downloaded them via `python -m sparc.download_data`, "
                    "or set $SPARC_PSP_PATH"
                ),
                color="FAIL",
            )
            return False


def api_test():
    cprint("Testing JSON API...", color="COMMENT")
    from sparc.api import SparcAPI

    try:
        api = SparcAPI()
    except Exception:
        return False
    version = api.sparc_version
    if version is None:
        cprint("Loaded API but no version date is provided.", color="WARNING")
    else:
        cprint(f"Loaded API version {version}", color="OKBLUE")
    return True


def command_test():
    cprint("Testing SPARC command...", color="COMMENT")
    from sparc.calculator import SPARC
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(directory=tmpdir)
        try:
            test_cmd = calc._make_command()
        except Exception:
            cprint(
                (
                    "No SPARC command found! \n"
                    "Please make sure you have sparc in your path, "
                    "or set up $ASE_SPARC_COMMAND variable!"
                ),
                color="FAIL",
            )
            return False
        cprint(f"The prefix for SPARC command is {test_cmd}", color="OKBLUE")
        return True


def calc_test():
    cprint("Running simple calculation...", color="COMMENT")
    from sparc.calculator import SPARC
    from ase.build import bulk
    import tempfile

    # 1x Al atoms with super bad calculation condition
    al = bulk("Al", cubic=False)

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.3, kpts=(1, 1, 1), tol_scf=1e-3, directory=tmpdir)
        try:
            al.calc = calc
            al.get_potential_energy()
        except Exception:
            cprint(
                ("Simple SPARC calculation failed!"),
                color="FAIL",
            )
            return False
        return True


def main():
    cprint(
        ("Performing a quick test on your " "SPARC and python API setup"),
        color=None,
    )
    results = {}
    results["Import"] = import_test()
    results["Pseudopotential"] = psp_test()
    results["JSON API"] = api_test()
    results["SPARC command"] = command_test()
    results["Calculation"] = False if results["SPARC command"] is False else calc_test()

    cprint(
        "\nSummary of test results",
        color="HEADER",
    )

    print_wiki = False
    for key, val in results.items():
        cprint(f"{key}:", bold=True, end="")
        if val is True:
            cprint("  PASS", color="OKGREEN")
        else:
            cprint("  FAIL", color="FAIL")
            print_wiki = True

    if print_wiki:
        cprint(
            "\nSome of the tests failed! Please refer to the following resources: \n"
            "1. SPARC's documentation: https://github.com/SPARC-X/SPARC/blob/master/doc/Manual.pdf \n"
            "2. Python API documentation: https://github.com/alchem0x2A/SPARC-X-API/blob/master/README.md\n",
            color="FAIL",
        )


if __name__ == "__main__":
    main()
