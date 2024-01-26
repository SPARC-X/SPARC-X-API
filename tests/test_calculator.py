import tempfile
from pathlib import Path

import pytest


def test_h_parameter():
    """Parameter h will be overwritten by any of FD_GRID, MESH_SPACING, ECUT"""
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "MESH_SPACING: 0.37" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, MESH_SPACING=0.4, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID:" not in filecontent
        assert "MESH_SPACING: 0.4" in filecontent  # Use SPARC mesh_spacing

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, ECUT=25, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "MESH_SPACING:" not in filecontent
        assert "FD_GRID:" not in filecontent
        assert "ECUT:" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, FD_GRID=[25, 25, 25], directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID: 25 25 25" in filecontent
        assert "MESH_GRID:" not in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(gpts=[25, 25, 25], directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID: 25 25 25" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.25, gpts=[25, 25, 25], directory=tmpdir)
        with pytest.raises(Exception):
            calc.write_input(atoms)

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(mesh_spacing=0.4, gpts=[25, 25, 25], directory=tmpdir)
        # gpts --> FD_GRID, this is excluded
        with pytest.raises(Exception):
            calc.write_input(atoms)

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(mesh_spacing=0.4, fd_grid=[25, 25, 25], directory=tmpdir)
        # gpts --> FD_GRID, this is excluded
        with pytest.raises(Exception):
            calc.write_input(atoms)


def test_xc_parameter():
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: GGA_PBE" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="pbe", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: GGA_PBE" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="lda", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: LDA_PZ" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="lda", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: LDA_PZ" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="rpbe", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: GGA_RPBE" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="pbesol", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: GGA_PBEsol" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="pbe0", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: PBE0" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="hf", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: HF" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="hse", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: HSE" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="hse03", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: HSE" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="vdw-df", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: vdWDF1" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="vdw-df2", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: vdWDF2" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(xc="scan", directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "EXCHANGE_CORRELATION: SCAN" in filecontent


def test_conflict_param():
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, FD_GRID=[25, 25, 25], ECUT=25)
        # FD_GRID and ECUT are conflict, but only detected during the writing
        with pytest.raises(Exception):
            calc.write_input(atoms)

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, FD_GRID=[25, 25, 25], MESH_SPACING=0.4)
        # FD_GRID and ECUT are conflict, but only detected during the writing
        with pytest.raises(Exception):
            calc.write_input(atoms)


def test_cell_param():
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, CELL=[10, 10, 10])
        # Cannot write the cell parameter if the atoms has already been provided
        with pytest.raises(Exception):
            calc.write_input(atoms)


def test_unknown_params():
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        with pytest.raises(Exception):
            calc = SPARC(h=0.2, directory=tmpdir, CUSTOM_KEY="[10, 10, 10]")


def test_label():
    from ase.build import bulk

    from sparc.calculator import SPARC

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, label="test_label")
        assert calc.label == "test_label"
        calc.write_input(atoms)
        assert (Path(tmpdir) / "test_label.inpt").exists()


def test_cache_results():
    # Test if the calculation results are cached (same structure)
    from pathlib import Path

    from ase.build import molecule

    from sparc.calculator import SPARC

    nh3 = molecule("NH3", cell=(8, 8, 8), pbc=True)
    nh3.rattle()

    dummy_calc = SPARC()
    try:
        cmd = dummy_calc._make_command()
    except EnvironmentError:
        print("Skip test since no sparc command found")
        return
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", print_forces=True, directory=tmpdir
        )
        nh3.calc = calc
        forces = nh3.get_forces()
        # make sure no more calculations are needed
        assert len(calc.check_state(nh3)) == 0
        energy = nh3.get_potential_energy()
        static_files = list(Path(tmpdir).glob("*.static*"))
        assert len(static_files) == 1


def test_sparc_version():
    from pathlib import Path

    from ase.build import molecule

    from sparc.calculator import SPARC

    h2 = molecule("H2", cell=[6, 6, 6], pbc=False)
    dummy_calc = SPARC()
    try:
        cmd = dummy_calc._make_command()
    except EnvironmentError:
        print("Skip test since no sparc command found")
        return
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", print_forces=True, directory=tmpdir
        )
        version = calc.detect_sparc_version()
        assert version is not None
        assert calc.directory == Path(tmpdir)
        assert calc.sparc_version is None
        calc1 = SPARC(
            h=0.3,
            kpts=(1, 1, 1),
            xc="pbe",
            print_forces=True,
            directory=tmpdir,
            check_version=True,
        )
        assert calc1.sparc_version is not None


def test_sparc_overwrite_files():
    """Check if keep_old_files option works"""
    from pathlib import Path

    from ase.build import molecule
    from ase.optimize import BFGS

    from sparc.calculator import SPARC

    h2 = molecule("H2", cell=[6, 6, 6], pbc=False)
    h2.center()
    dummy_calc = SPARC()
    try:
        cmd = dummy_calc._make_command()
    except EnvironmentError:
        print("Skip test since no sparc command found")
        return

    # Default, keep_old_file=False, only 1 static file
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", directory=tmpdir, keep_old_files=False
        )
        atoms = h2.copy()
        atoms.calc = calc
        opt = BFGS(atoms)
        opt.run(fmax=2.0)
        assert len(list(tmpdir.glob("*.static*"))) == 1
        assert len(list(tmpdir.glob("*.out*"))) == 1

    # Default, keep_old_file=True, multiple files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", directory=tmpdir, keep_old_files=True
        )
        atoms = h2.copy()
        atoms.calc = calc
        opt = BFGS(atoms)
        opt.run(fmax=2.0)
        # On ase 3.22 with default BFGS parameters, it taks 6 steps
        assert len(list(tmpdir.glob("*.static*"))) > 3
        assert len(list(tmpdir.glob("*.out*"))) > 3


def test_incompat_atoms():
    """SPARC calculator should capture incorrect atoms positions before writing
    the inputs
    """
    from ase.build import bulk, molecule

    from sparc.calculator import SPARC

    calc = SPARC()

    # CASE 1: Zero-length cell
    atoms = molecule("CH4")
    with pytest.raises(ValueError):
        calc.calculate(atoms)

    # CASE 2: PBC wrapp back should be valid
    atoms = molecule("CH4", cell=[6, 6, 6], pbc=True)  # some H atoms are below 0
    assert atoms.positions[:, 0].min() < 0
    assert calc.check_input_atoms(atoms) is None

    # CASE 3: Dirichlet BC and atoms outside domain
    atoms = molecule("CH4", cell=[6, 6, 6], pbc=False)  # some H atoms are below 0
    assert atoms.positions[:, 0].min() < 0
    with pytest.raises(ValueError):
        calc.check_input_atoms(atoms)

    atoms.center()  # Centered atoms should be allowed
    assert calc.check_input_atoms(atoms) is None

    atoms.positions[:, 2] += 4  # Move up the z-direction until atoms are out-of-domain
    with pytest.raises(ValueError):
        calc.check_input_atoms(atoms)

    # CASE 4: Non-orthogonal cell
    atoms = bulk("Al", cubic=False)
    assert atoms.cell.angles()[0] != 90.0
    print(atoms.pbc)
    assert calc.check_input_atoms(atoms) is None  # Primitive cell is ok for pbc=True
    atoms.pbc = False  # Unphysical structure, just for checking
    with pytest.raises(ValueError):
        calc.check_input_atoms(atoms)


def test_calc_param_set():
    """Test if the set method works"""
    from sparc.calculator import SPARC

    # CASE 1: default params
    calc = SPARC()
    assert calc.label == "SPARC"  # Non-restart, use default label
    assert calc.valid_params == {}
    assert set(calc.special_params.keys()) == set(["xc", "kpts", "h"])

    calc.set(h=0.22)
    assert "h" in calc.special_params

    calc.set(gpts=[50, 50, 50])
    assert "h" not in calc.special_params
    assert "gpts" in calc.special_params

    calc.set(calc_stress=False)
    assert "CALC_STRESS" in calc.valid_params
    assert "CALC_STRESS" in calc.parameters

    with pytest.raises(KeyError):
        # Invalid parameter
        calc.set(calc_stres=False)

    with pytest.raises(ValueError):
        calc.set(calc_stress="ok")

    calc.set(directory="test")
    assert calc.directory.name == "test"
    assert "directory" not in calc.parameters

    # CASE 2: change order of parameters
    calc = SPARC(xc="pbe", h=0.22, directory="test", log="test.log")
    assert "directory" not in calc.parameters
    assert "log" not in calc.parameters
    assert str(calc.log) == "test/test.log"


def test_sparc_param_change():
    """Change of either inpt_state or system_state would case
    calculator to reload
    """
    from pathlib import Path

    from ase.build import molecule
    from ase.optimize import BFGS

    from sparc.calculator import SPARC

    h2 = molecule("H2", cell=[6, 6, 6], pbc=False)
    h2.center()
    dummy_calc = SPARC()
    try:
        cmd = dummy_calc._make_command()
    except EnvironmentError:
        print("Skip test since no sparc command found")
        return

    # Default, keep_old_file=False, only 1 static file
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", directory=tmpdir, keep_old_files=False
        )
        atoms = h2.copy()
        atoms.calc = calc
        opt = BFGS(atoms)
        opt.run(fmax=2.0)
        assert len(list(tmpdir.glob("*.static*"))) == 1
        assert len(list(tmpdir.glob("*.out*"))) == 1

    # Default, keep_old_file=True, multiple files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(
            h=0.3, kpts=(1, 1, 1), xc="pbe", directory=tmpdir, keep_old_files=True
        )
        atoms = h2.copy()
        atoms.calc = calc
        opt = BFGS(atoms)
        opt.run(fmax=2.0)
        # On ase 3.22 with default BFGS parameters, it taks 6 steps
        assert len(list(tmpdir.glob("*.static*"))) > 3
        assert len(list(tmpdir.glob("*.out*"))) > 3
