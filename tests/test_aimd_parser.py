import os
from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr, Hartree, eV, kB

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_aimd_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.aimd import _read_aimd

    data_dict = _read_aimd(
        test_output_dir
        / "TiO2_orthogonal_quick_md.sparc"
        / "TiO2_orthogonal_quick_md.aimd"
    )
    assert "aimd" in data_dict
    md_steps = data_dict["aimd"]
    assert len(md_steps) == 5
    for i, step in enumerate(md_steps):
        assert i == step.get("step", -1)
        for key in [
            "positions",
            "velocities",
            "forces",
            "electron temp",
            "ion temp",
            "total energy per atom",
            "kinetic energy per atom",
            "kinetic energy (ideal gas) per atom",
            "free energy per atom",
            "stress",
            "pressure",
        ]:
            assert key in step
        assert step["positions"].shape == (6, 3)
        assert step["forces"].shape == (6, 3)
        assert step["velocities"].shape == (6, 3)
        # TODO: may subject to changes
        assert step["stress"].shape == (3, 3)

    # A simple test to see if we're using the correct types
    T0 = md_steps[0]["electron temp"]
    ek_ig0 = md_steps[0]["kinetic energy (ideal gas) per atom"]
    assert np.isclose(T0, 800)
    assert np.isclose(T0, ek_ig0 * eV / kB / (1.5), 1.0e-2)
