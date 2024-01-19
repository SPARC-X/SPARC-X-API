import tempfile
from pathlib import Path

import pytest

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"
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

def test_socket_read_same_cell():
    """Test parsing multi-image static files from socket"""
    from sparc.io import read_sparc
    from ase.io import read
    import numpy as np
    bundle = test_output_dir / "Al_socket_bfgs.sparc"
    # References from standard socket output
    ref_images = read(bundle / "sparc-socket.traj", ":")
    cell0 = ref_images[0].cell
    parsed_images = read_sparc(bundle, ":")
    assert len(ref_images) == len(parsed_images)
    for r_img, p_img in zip(ref_images, parsed_images):
        assert np.isclose(r_img.positions, p_img.positions, 1.e-4).all()
        assert np.isclose(r_img.cell, p_img.cell, 1.e-4).all()
        assert np.isclose(cell0, p_img.cell, 1.e-4).all()
        assert np.isclose(r_img.get_potential_energy(), p_img.get_potential_energy(), 1.e-4)
        assert np.isclose(r_img.get_forces(), p_img.get_forces(), 1.e-4).all()


def test_socket_read_diff_cell():
    """Test parsing multi-image static files from socket. The volumes of the cells change"""
    from sparc.io import read_sparc
    import numpy as np
    bundle = test_output_dir / "Al_socket_volchange.sparc"
    parsed_images = read_sparc(bundle, ":")
    assert len(parsed_images) == 10
    cell0 = parsed_images[0].cell

    for i, p_img in enumerate(parsed_images):
        bundle_ref = bundle / "single-points" / f"sp_image{i:02d}"
        r_img = read_sparc(bundle_ref)
        p_img.wrap()
        assert np.isclose(r_img.positions, p_img.positions, 1.e-4).all()
        assert np.isclose(r_img.cell, p_img.cell, 1.e-4).all()
        assert np.isclose(cell0, p_img.cell, 1.e-4).all()
        assert np.isclose(r_img.get_potential_energy(), p_img.get_potential_energy(), 1.e-4)
        assert np.isclose(r_img.get_forces(), p_img.get_forces(), 1.e-4).all()
