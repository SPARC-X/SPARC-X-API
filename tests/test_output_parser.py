import pytest
import numpy as np
from pathlib import Path
import os
from ase.units import Bohr, Hartree, kB, eV

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_output_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.out import _read_out

    data_dict = _read_out(
        test_output_dir / "AlSi_primitive_quick_relax.sparc" /
        "AlSi_primitive_quick_relax.out"
    )
    assert "out" in data_dict
    out_dict = data_dict["out"]
    assert out_dict["sparc_version"] == "2023.02.03"
    # Currently the fields in run_info are strings
    assert int(out_dict["run_info"]["number of processors"]) == 48
    assert out_dict["parameters"]["NP_SPIN_PARAL"] == 1
    assert out_dict["parameters"]["LATVEC"].shape == (3, 3)
    ionic_steps = out_dict["ionic_steps"]
    assert len(ionic_steps) == 5
    for step in ionic_steps:
        assert "convergence" in step
        assert "total free energy" in step
