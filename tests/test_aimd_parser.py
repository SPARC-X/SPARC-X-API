import pytest
import numpy as np
from pathlib import Path
import os
from ase.units import Bohr, Hartree

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_aimd_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.aimd import _read_aimd

    data_dict = _read_aimd(test_output_dir / "TiO2_orthogonal_quick_md.sparc" / "TiO2_orthogonal_quick_md.aimd")
    assert "aimd" in data_dict
    md_steps = data_dict["aimd"]
    for i, step in enumerate(md_steps):
        assert i == step.get("step", -1)
        for key in ["positions", "velocities", "forces", "electron temp",
                    "total energy per atom", "kinetic energy per atom",
                    "free energy per atom", "stress", "pressure"]:
            assert key in step
