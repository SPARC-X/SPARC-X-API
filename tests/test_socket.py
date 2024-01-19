import tempfile
from pathlib import Path

import pytest

curdir = Path(__file__).parent
repo_dir = curdir.parent


socket_stdout = b"""
USAGE:
    mpirun -np <nproc> {SPARCROOT}/lib/sparc -name <filename>

    {SPARCROOT} is the location of the SPARC folder

REQUIRED ARGUMENT:
    -name <filename>
           The filename shared by .inpt file and .ion
           file (without extension)

OPTIONS: 
    -h, --help
           Display help (from command line).
    -n <number of Nodes>
    -c <number of CPUs per node>
    -a <number of Accelerators (e.g., GPUs) per node>
    -socket <socket> 
            <socket> can be either  <host>:<port> or <unix_socket>:UNIX.
            Note: socket (driver) mode is an experimental feature.

EXAMPLE:

    mpirun -np 8 {SPARCROOT}/lib/sparc -name test

    The example command runs sparc with 8 cores, with input file named
    test.inpt, and ion file named test.ion.

NOTE: 
    This is a short description of the usage of SPARC. For a detailed 
    discription, refer to the manual online at

        https://github.com/SPARC-X/SPARC/tree/master/doc
"""

non_socket_stdout = b"""
USAGE:
    mpirun -np <nproc> {SPARCROOT}/lib/sparc -name <filename>

    {SPARCROOT} is the location of the SPARC folder

REQUIRED ARGUMENT:
    -name <filename>
           The filename shared by .inpt file and .ion
           file (without extension)

OPTIONS: 
    -h, --help
           Display help (from command line).
    -n <number of Nodes>
    -c <number of CPUs per node>
    -a <number of Accelerators (e.g., GPUs) per node>

EXAMPLE:

    mpirun -np 8 {SPARCROOT}/lib/sparc -name test

    The example command runs sparc with 8 cores, with input file named
    test.inpt, and ion file named test.ion.

NOTE: 
    This is a short description of the usage of SPARC. For a detailed 
    discription, refer to the manual online at

        https://github.com/SPARC-X/SPARC/tree/master/doc 
"""


def test_socket_compat(monkeypatch):
    """Test if SPARC calculator can recognize socket version binary"""
    from ase.build import bulk

    from sparc.calculator import SPARC

    def mock_nonsocket_run(*args, **kwargs):
        MockCompletedProcess = type("CompletedProcess", (object,), {})
        mock_cp = MockCompletedProcess()
        mock_cp.stdout = non_socket_stdout
        return mock_cp

    def mock_socket_run(*args, **kwargs):
        MockCompletedProcess = type("CompletedProcess", (object,), {})
        mock_cp = MockCompletedProcess()
        mock_cp.stdout = socket_stdout
        return mock_cp

    calc = SPARC(command="sparc")
    # Mock a non-socket sparc binary stdout
    monkeypatch.setattr("subprocess.run", mock_nonsocket_run)
    assert calc.detect_socket_compatibility() is False

    # Mock a socket sparc binary stdout
    monkeypatch.setattr("subprocess.run", mock_socket_run)
    assert calc.detect_socket_compatibility() is True
