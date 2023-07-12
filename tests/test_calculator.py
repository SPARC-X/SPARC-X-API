import pytest
from pathlib import Path
import tempfile


def test_h_parameter():
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
    
