import os
from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr, Hartree

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_geopt_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.geopt import _read_geopt

    data_dict = _read_geopt(
        test_output_dir
        / "AlSi_primitive_quick_relax.sparc"
        / "AlSi_primitive_quick_relax.geopt"
    )
    assert "geopt" in data_dict
    geopt_steps = data_dict["geopt"]
    for i, step in enumerate(geopt_steps):
        assert i == step.get("step", -1)
        assert "positions" in step
        assert "forces" in step
        assert "energy" in step
        assert "cell" in step
        assert "volume" in step
        assert "latvec" in step
        assert "stress" in step
        assert "ase_cell" in step

        # Value assertions
        ase_cell = step["ase_cell"]
        vol_ase = np.linalg.det(ase_cell)
        assert np.isclose(vol_ase, step["volume"])

    max_final_f = np.max(np.abs(step["forces"]))
    assert max_final_f < 1.0e-3 * Hartree / Bohr


def test_geopt_parser_relax2():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.geopt import _read_geopt

    data_dict = _read_geopt(
        test_output_dir / "Si8_cell_geopt_relax2.sparc" / "Si8_cell_geopt.geopt"
    )
    geopt_steps = data_dict["geopt"]
    for i, step in enumerate(geopt_steps):
        assert i == step.get("step", -1)
        # RELAX=2 no position information
        assert "positions" not in step
        assert "stress" in step
        assert "cell" in step
        assert "latvec" in step


def test_geopt_low_dim_stress():
    from sparc.sparc_parsers.geopt import _read_geopt

    data_dict = _read_geopt(
        test_output_dir / "Alloy_geopt_ppd_bc.sparc" / "SPARC.geopt"
    )
    geopt_steps = data_dict["geopt"]
    for i, step in enumerate(geopt_steps):
        assert i == step.get("step", -1)
        assert "stress_2d" in step
        assert "stress" not in step
        assert "cell" in step
