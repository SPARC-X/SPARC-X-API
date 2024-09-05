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
    import numpy as np
    from ase.io import read

    from sparc.io import read_sparc

    bundle = test_output_dir / "Al_socket_bfgs.sparc"
    # References from standard socket output
    ref_images = read(bundle / "sparc-socket.traj", ":")
    cell0 = ref_images[0].cell
    parsed_images = read_sparc(bundle, ":")
    assert len(ref_images) == len(parsed_images)
    for r_img, p_img in zip(ref_images, parsed_images):
        assert np.isclose(
            r_img.get_potential_energy(), p_img.get_potential_energy(), 1.0e-4
        )
        assert np.isclose(r_img.get_forces(), p_img.get_forces(), 1.0e-4).all()
        assert np.isclose(r_img.positions, p_img.positions, 1.0e-4).all()
        assert np.isclose(r_img.cell, p_img.cell, 1.0e-4).all()
        assert np.isclose(cell0, p_img.cell, 1.0e-4).all()


def test_socket_read_diff_cell():
    """Test parsing multi-image static files from socket. The volumes of the cells change"""
    import numpy as np

    from sparc.io import read_sparc

    bundle = test_output_dir / "Al_socket_volchange.sparc"
    parsed_images = read_sparc(bundle, ":")
    assert len(parsed_images) == 10
    cell0 = parsed_images[0].cell

    ratios = [1.0, 1.01, 1.02, 1.03, 1.04, 1.0, 0.99, 0.98, 0.97, 0.96]
    for i, p_img in enumerate(parsed_images):
        print(i)
        bundle_ref = bundle / "single-points" / f"sp_image{i:02d}"
        r_img = read_sparc(bundle_ref)
        # For systems with cell change, force error may be larger. This may be due to mesh resizing
        assert np.isclose(
            r_img.get_potential_energy(), p_img.get_potential_energy(), rtol=1.0e-3
        )
        assert np.isclose(r_img.get_forces(), p_img.get_forces(), atol=6.0e-2).all()
        # Wrap may be necessary for p_img
        p_img.wrap()
        assert np.isclose(r_img.positions, p_img.positions, 1.0e-4).all()
        assert np.isclose(r_img.cell, p_img.cell, 1.0e-4).all()
        assert np.isclose(cell0 * ratios[i], p_img.cell, 1.0e-4).all()


def test_socket_incomplete():
    """Test an incomplete socket calculation output
    The test sample is problematic, but for the purpose of demonstrating
    reading incomplete files should be enough
    """
    from sparc.io import read_sparc

    bundle = test_output_dir / "H2O_socket_incomplete.sparc"
    parsed_images = read_sparc(bundle, ":")
    assert len(parsed_images) == 4
    for i in range(3):
        assert parsed_images[i].get_potential_energy() is not None
        assert parsed_images[i].get_forces() is not None
    with pytest.raises(Exception):
        parsed_images[3].get_potential_energy()
    with pytest.raises(Exception):
        parsed_images[3].get_forces()


def test_socket_param_calculator():
    """Test if socket settings are correct"""
    from sparc.calculator import SPARC

    # Default use_socket is False
    calc = SPARC()
    assert calc.use_socket is False

    # CASE 1: use alias to switch on socket mode
    calc = SPARC(use_socket=True)
    assert calc.use_socket
    assert calc.socket_params["use_socket"]
    assert calc.socket_params["host"] == "localhost"
    assert calc.socket_params["port"] == -1
    assert calc.socket_params["server_only"] is False

    # CASE 2: just use socket params. socket params
    # have higher priority
    calc = SPARC(socket_params=dict(use_socket=True))
    assert calc.use_socket

    # CASE 3: partial change socket params
    calc = SPARC(socket_params=dict(host="server"))
    assert calc.use_socket is False

    # CASE 4: partial change socket params
    calc = SPARC(use_socket=True, socket_params=dict(host="server"))
    assert calc.use_socket is True
    assert calc.socket_params["host"] == "server"


def test_atoms_system_changes(monkeypatch):
    """Test if change of atoms will cause socket redraw"""
    from ase.atoms import Atom
    from ase.build import molecule

    from sparc.calculator import SPARC, all_changes

    def mock_calculate(
        self, atoms=None, properties=["energy"], system_changes=all_changes
    ):
        # Mock implementation
        self.results["energy"] = 0.0
        self.results["forces"] = [[0, 0, 0]] * len(atoms)
        self.atoms = atoms.copy()
        self.atoms.arrays["initial_magmoms"] = [
            0,
        ] * len(atoms)

    monkeypatch.setattr(SPARC, "_calculate_with_socket", mock_calculate)

    calc = SPARC(use_socket=True)
    # calc._calculate_with_socket = mock_calculate

    # Case 1: change of positions
    h2o = molecule("H2O", cell=[6, 6, 6], pbc=True)
    h2o.calc = calc
    h2o.get_potential_energy()
    h2o_new = h2o.copy()
    h2o_new.rattle()
    # We can pop out "system_state" because these are fake calculations
    changes = set(calc.check_state(h2o_new))
    changes.discard("system_state")
    # h2o_new.calc = calc
    # h2o_new.get_potential_energy()
    assert set(changes) == set(["positions"])

    # Case 2: change of positions and cell
    h2o_new.set_cell([8, 8, 8], scale_atoms=True)
    changes = set(calc.check_state(h2o_new))
    changes.discard("system_state")
    assert set(changes) == set(["positions", "cell"])

    # Case 3: change of pbc
    h2o_new = h2o.copy()
    h2o_new.pbc = [False, True, False]  # arbitrary change
    changes = set(calc.check_state(h2o_new))
    changes.discard("system_state")
    assert set(changes) == set(["pbc"])

    # Case 4: change of numbers
    h2o_new = h2o.copy()
    h2o_new += Atom("H", [5, 5, 5])
    changes = set(calc.check_state(h2o_new))
    changes.discard("system_state")
    assert "numbers" in changes


def test_atoms_system_changes_real():
    """Real system changes test"""
    from ase.build import molecule

    from sparc.calculator import SPARC, all_changes

    h2o = molecule("H2O", cell=[6, 6, 6], pbc=True)
    h2o.center()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(use_socket=True, directory="test-change", h=0.25)
        try:
            can_run = calc.detect_socket_compatibility()
        except Exception:
            can_run = False
        if not can_run:
            print("Cannot start a sparc calculation with socket. Skip test")
            pytest.skip()

        assert calc.in_socket is not None
        infile0 = calc.in_socket_filename
        assert infile0.startswith("/tmp/ipi_sparc_")
        calc.close()
        assert calc.in_socket is None
        h2o.calc = calc
        h2o.get_potential_energy()
        infile1 = calc.in_socket_filename
        assert infile0 != infile1
        assert calc.pid is not None
        old_pid = calc.pid
        # Case 1: positions change, same socket process
        h2o.rattle()
        h2o.get_potential_energy()
        assert calc.pid == old_pid

        # Case 2: cell and position change, same socket
        h2o.set_cell([6.05, 6.05, 6.05], scale_atoms=True)
        h2o.get_potential_energy()
        assert calc.pid == old_pid

        # Case 3: pbc change, restart
        h2o.pbc = [False, False, False]
        h2o.get_potential_energy()
        assert calc.pid != old_pid

        # Case 4: change of calc parameters, this will restart
        calc.set(h=0.27)
        h2o.get_potential_energy()
        assert calc.pid != old_pid

        # Case 4: change of system state, will restart
        calc.set(label="SPARC-1")
        h2o.get_potential_energy()
        assert calc.pid != old_pid

        calc.close()
