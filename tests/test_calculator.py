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
