import pytest
from pathlib import Path
import tempfile


def test_h_parameter():
    """Parameter h will be overwritten by any of FD_GRID, MESH_SPACING, ECUT"""
    from sparc.calculator import SPARC
    from ase.build import bulk

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID: 21 21 21" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, MESH_SPACING=0.4, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID:" not in filecontent
        assert "MESH_SPACING: 0.4" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, ECUT=25, directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID:" not in filecontent
        assert "ECUT:" in filecontent

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, FD_GRID=[25, 25, 25], directory=tmpdir)
        calc.write_input(atoms)
        filecontent = open(Path(tmpdir) / "SPARC.inpt", "r").read()
        assert "FD_GRID: 25 25 25" in filecontent


def test_conflict_param():
    from sparc.calculator import SPARC
    from ase.build import bulk

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, FD_GRID=[25, 25, 25], ECUT=25)
        # FD_GRID and ECUT are conflict, but only detected during the writing
        with pytest.raises(Exception):
            calc.write_input(atoms)

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(
            h=0.2, directory=tmpdir, FD_GRID=[25, 25, 25], MESH_SPACING=0.4
        )
        # FD_GRID and ECUT are conflict, but only detected during the writing
        with pytest.raises(Exception):
            calc.write_input(atoms)


def test_cell_param():
    from sparc.calculator import SPARC
    from ase.build import bulk

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, CELL=[10, 10, 10])
        # Cannot write the cell parameter if the atoms has already been provided
        with pytest.raises(Exception):
            calc.write_input(atoms)


def test_unknown_params():
    from sparc.calculator import SPARC
    from ase.build import bulk

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        with pytest.raises(Exception):
            calc = SPARC(h=0.2, directory=tmpdir, CUSTOM_KEY="[10, 10, 10]")


def test_label():
    from sparc.calculator import SPARC
    from ase.build import bulk

    atoms = bulk("Al", cubic=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = SPARC(h=0.2, directory=tmpdir, label="test_label")
        assert calc.label == "test_label"
        calc.write_input(atoms)
        assert (Path(tmpdir) / "test_label.inpt").exists()


def test_cache_results():
    # Test if the calculation results are cached (same structure)
    from sparc.calculator import SPARC
    from ase.build import molecule
    from pathlib import Path

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
