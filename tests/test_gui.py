import os
from pathlib import Path
from subprocess import run

import numpy as np
import pytest
from ase.units import Bohr, Hartree, eV, kB

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_cli():
    """Simply call the sparc-ase methods"""
    proc = run(["sparc-ase"])
    assert proc.returncode == 0


def test_info():
    """Call the sparc-ase info on sparc file. Be aware of the API changes in 3.22->3.23"""
    import ase
    from packaging import version

    bundle = test_output_dir / "Cu_FCC.sparc"
    if version.parse(ase.__version__) < version.parse("3.23"):
        proc = run(["sparc-ase", "info", f"{bundle}"], capture_output=True)
    else:
        proc = run(["sparc-ase", "info", "--files", f"{bundle}"], capture_output=True)
    assert proc.returncode == 0
    assert "SPARC" in proc.stdout.decode("utf8")


def test_gui_single():
    """Call the sparc-ase gui on a static calculation"""
    bundle = test_output_dir / "Cu_FCC.sparc"
    proc = run(["sparc-ase", "gui", "-t", f"{bundle}"], capture_output=True)
    assert proc.returncode == 0


def test_gui_geopt():
    """Call the sparc-ase info on a geopt calculation"""
    bundle = test_output_dir / "AlSi_primitive_quick_relax.sparc"
    proc = run(["sparc-ase", "gui", "-t", f"{bundle}"], capture_output=True)
    assert proc.returncode == 0


def test_gui_md():
    """Call the sparc-ase info on an AIMD calculation"""
    bundle = test_output_dir / "TiO2_orthogonal_quick_md.sparc"
    proc = run(["sparc-ase", "gui", "-t", f"{bundle}"], capture_output=True)
    assert proc.returncode == 0
