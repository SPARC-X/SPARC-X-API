"""A simple test module for sparc python api
Usage:
python -m sparc.quicktest
"""
from pathlib import Path

from ase.data import chemical_symbols

from .utils import cprint


class BaseTest(object):
    """Base class for all tests providing functionalities

    Each child class will implement its own `run_test` method to
    update the `result`, `error_handling` and `info` fields.

    If you wish to include a simple error handling message for each
    child class, add a line starting `Error handling` follows by the
    helper message at the end of the docstring
    """

    def __init__(self):
        self.result = None
        self.error_msg = ""
        self.error_handling = ""
        self.info = {}

    @property
    @classmethod
    def dislay_name(cls):
        return cls.__name__

    def display_docstring(self):
        """Convert the class's docstring to error handling"""
        doc = self.__class__.__doc__
        error_handling_lines = []
        begin_record = False
        indent = 0  # indentation for the "Error handling" line
        if doc:
            for line in doc.splitlines():
                if line.lstrip().startswith("Error handling"):
                    if begin_record is True:
                        msg = (
                            "There are multiple Error handlings "
                            "in the docstring of "
                            f"{self.__class__.__name__}."
                        )
                        raise ValueError(msg)
                    begin_record = True
                    indent = len(line) - len(line.lstrip())
                elif begin_record is True:
                    current_indent = len(line) - len(line.lstrip())
                    line = line.strip()
                    if len(line) > 0:  # Only add non-empty lines
                        # Compensate for the extra indentation
                        # if current_indent > indent
                        spaces = max(0, current_indent - indent) * " "
                        error_handling_lines.append(spaces + line)
                else:
                    pass
        else:
            pass
        error_handling_string = "\n".join(error_handling_lines)
        return error_handling_string

    def make_test(self):
        """Each class should implement ways to update `result` and `info`"""
        raise NotImplementedError

    def run_test(self):
        """Run test and update result etc.
        If result is False, update the error handling message
        """
        try:
            self.make_test()
        except Exception as e:
            self.result = False
            self.error_msg = str(e)

        if self.result is None:
            raise ValueError(
                "Test result is not updated for " f"{self.__class__.__name__} !"
            )
        if self.result is False:
            self.error_handling = self.display_docstring()
        return


class ImportTest(BaseTest):
    """Check if external io format `sparc` can be registered in ASE

    Error handling:
    - Make sure SPARC-X-API is installed via conda / pip / setuptools
    - If you wish to work on SPARC-X-API source code, use `pip install -e`
      instead of setting up $PYTHON_PATH
    """

    display_name = "Import"

    def make_test(self):
        cprint("Testing import...", color="COMMENT")
        from ase.io.formats import ioformats

        self.result = "sparc" in ioformats.keys()
        if self.result is False:
            self.error_msg = (
                "Cannot find `sparc` as a valid " "external ioformat for ASE."
            )
        return


class PspTest(BaseTest):
    """Check at least one directory of Pseudopotential files exist
    info[`psp_dir`] contains the first psp dir found on system
    # TODO: check if all psp files can be located
    #TODO: update to the ASE 3.23 config method

    Error handling:
    - Default version of psp files can be downloaded by
      `python -m sparc.download_data`
    - Alternatively, specify the variable $SPARC_PSP_PATH
      to the custom pseudopotential files
    """

    display_name = "Pseudopotential"

    def make_test(self):
        cprint("Testing pseudo potential path...", color="COMMENT")
        import tempfile

        from .io import SparcBundle
        from .sparc_parsers.pseudopotential import find_pseudo_path

        with tempfile.TemporaryDirectory() as tmpdir:
            sb = SparcBundle(directory=tmpdir)
            psp_dir = sb.psp_dir

        if psp_dir is not None:
            psp_dir = Path(psp_dir)
            self.info["psp_dir"] = f"{psp_dir.resolve()}"
            if not psp_dir.is_dir():
                self.result = False
                self.error_msg = (
                    "Pseudopotential files path " f"{psp_dir.resolve()} does not exist."
                )
            else:
                missing_elements = []
                # Default psp file are 1-57 + 72-83
                spms_elements = chemical_symbols[1:58] + chemical_symbols[72:84]
                for element in spms_elements:
                    try:
                        find_pseudo_path(element, psp_dir)
                    except Exception:
                        missing_elements.append(element)
                if len(missing_elements) == 0:
                    self.result = True
                else:
                    self.result = False
                    self.error_msg = (
                        "Pseudopotential files for "
                        f"{len(missing_elements)} elements are "
                        "missing or incompatible: \n"
                        f"{missing_elements}"
                    )
        else:
            self.info["psp_dir"] = "None"
            self.result = False
            self.error_msg = (
                "Pseudopotential file path not defined and/or "
                "default psp files are incomplete."
            )
        return


class ApiTest(BaseTest):
    """Check if the API can be loaded, and store the Schema version.

    # TODO: consider change to schema instead of api
    # TODO: allow config to change json file path
    Error handling:
    - Check if default JSON schema exists in
      `<sparc-x-api-root>/sparc_json_api/parameters.json`
    - Use $SPARC_DOC_PATH to specify the raw LaTeX files
    """

    display_name = "JSON API"

    def make_test(self):
        from .utils import locate_api

        try:
            api = locate_api()
            version = api.sparc_version
            self.result = True
            self.info["api_version"] = version
            self.info["api_source"] = api.source
        except Exception as e:
            self.result = False
            self.info["api_version"] = "NaN"
            self.info["api_source"] = "not found"
            self.error_msg = (
                "Error when locating a JSON schema or "
                f"LaTeX source files for SPARC. Error is {e}"
            )
        return


class CommandTest(BaseTest):
    """Check validity of command to run SPARC calculation. This test
    also checks sparc version and socket compatibility

    # TODO: check ase 3.23 config with separate binary
    Error handling:
    - The command prefix to run SPARC calculation should look like
      `<mpi instructions> <sparc binary>`
    - Use $ASE_SPARC_COMMAND to set the command string
    - Check HPC resources and compatibility (e.g. `srun` on a login node)
    """

    display_name = "SPARC Command"

    def make_test(self):
        import tempfile

        from sparc.calculator import SPARC

        self.info["command"] = ""
        self.info["sparc_version"] = ""

        with tempfile.TemporaryDirectory() as tmpdir:
            calc = SPARC(directory=tmpdir)
            # Step 1: validity of sparc command
            try:
                test_cmd = calc._make_command()
                self.result = True
                self.info["command"] = test_cmd
            except Exception as e:
                self.result = False
                self.info["command"] = "not found"
                self.error_msg = f"Error setting SPARC command:\n{e}"

            # Step 2: check SPARC binary version
            try:
                sparc_version = calc.detect_sparc_version()
                # Version may be None if failed to retrieve
                if sparc_version:
                    self.result = self.result & True
                    self.info["sparc_version"] = sparc_version
                else:
                    self.result = False
                    self.info["sparc_version"] = "NaN"
                    self.error_msg += "\n" if len(self.error_msg) > 0 else ""
                    self.error_msg += "Error detecting SPARC version"
            except Exception as e:
                self.result = False
                self.info["sparc_version"] = "NaN"
                self.error_msg += "\n" if len(self.error_msg) > 0 else ""
                self.error_msg += f"\nError detecting SPARC version:\n{e}"
        return


class FileIOCalcTest(BaseTest):
    """Run a simple calculation in File IO mode.

    # TODO: check ase 3.23 config
    Error handling:
    - Check if settings for pseudopotential files are correct
    - Check if SPARC binary exists and functional
    - Check if specific HPC requirements are met:
      (module files, libraries, parallel settings, resources)
    """

    display_name = "Calculation (File I/O)"

    def make_test(self):
        import tempfile

        from ase.build import bulk

        from sparc.calculator import SPARC

        # 1x Al atoms with super bad calculation condition
        al = bulk("Al", cubic=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            calc = SPARC(h=0.3, kpts=(1, 1, 1), tol_scf=1e-3, directory=tmpdir)
            try:
                al.calc = calc
                al.get_potential_energy()
                self.result = True
            except Exception as e:
                self.result = False
                self.error_msg = "Simple calculation in file I/O mode failed: \n" f"{e}"
        return


class SocketCalcTest(BaseTest):
    """Run a simple calculation in Socket mode (UNIX socket).

    # TODO: check ase 3.23 config
    Error handling:
    - The same as error handling in file I/O calculation test
    - Check if SPARC binary supports socket
    """

    display_name = "Calculation (UNIX socket)"

    def make_test(self):
        import tempfile

        from ase.build import bulk

        from sparc.calculator import SPARC

        # Check SPARC binary socket compatibility
        with tempfile.TemporaryDirectory() as tmpdir:
            calc = SPARC(directory=tmpdir)
            try:
                sparc_compat = calc.detect_socket_compatibility()
                self.info["sparc_socket_compatibility"] = sparc_compat
            except Exception:
                self.info["sparc_socket_compatibility"] = False

        # 1x Al atoms with super bad calculation condition
        al = bulk("Al", cubic=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            calc = SPARC(
                h=0.3, kpts=(1, 1, 1), tol_scf=1e-3, use_socket=True, directory=tmpdir
            )
            try:
                al.calc = calc
                al.get_potential_energy()
                self.result = True
            except Exception as e:
                self.result = False
                self.error_msg = (
                    "Simple calculation in socket mode (UNIX socket) failed: \n" f"{e}"
                )
        return


def main():
    cprint(
        ("Performing a quick test on your " "SPARC and python API setup"),
        color=None,
    )

    test_classes = [
        ImportTest(),
        PspTest(),
        ApiTest(),
        CommandTest(),
        FileIOCalcTest(),
        SocketCalcTest(),
    ]

    system_info = {}
    for test in test_classes:
        test.run_test()
        system_info.update(test.info)

    # Header section
    print("-" * 80)
    cprint(
        "Summary",
        bold=True,
        color="HEADER",
    )
    print("-" * 80)
    cprint("Configuration", bold=True, color="HEADER")
    for key, val in system_info.items():
        print(f"{key}: {val}")

    print("-" * 80)
    # Body section
    cprint("Tests", bold=True, color="HEADER")

    print_wiki = False
    for test in test_classes:
        cprint(f"{test.display_name}:", bold=True, end="")
        if test.result is True:
            cprint("  PASS", color="OKGREEN")
        else:
            cprint("  FAIL", color="FAIL")
            print_wiki = True

    print("-" * 80)
    # Error information section
    has_print_error_header = False
    for test in test_classes:
        if (test.result is False) and (test.error_handling):
            if has_print_error_header is False:
                cprint(
                    ("Some tests failed! " "Please check the following information.\n"),
                    color="FAIL",
                )
                has_print_error_header = True
            cprint(f"{test.display_name}:", bold=True)
            cprint(f"{test.error_msg}", color="FAIL")
            print(test.error_handling)
            print("\n")

    if print_wiki:
        print("-" * 80)
        cprint(
            "Please check additional information from:\n"
            "1. SPARC's documentation: https://github.com/SPARC-X/SPARC/blob/master/doc/Manual.pdf \n"
            "2. Python API documentation: https://github.com/alchem0x2A/SPARC-X-API/blob/master/README.md\n",
            color=None,
        )


if __name__ == "__main__":
    main()
