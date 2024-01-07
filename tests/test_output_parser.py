import os
from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr, Hartree, eV, kB

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_output_date_parser():
    from sparc.sparc_parsers.out import _read_sparc_version

    header1 = """***************************************************************************
*                       SPARC (version Feb 03, 2023)                      *
*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *
*           Distributed under GNU General Public License 3 (GPL)          *
*                   Start time: Sun Feb  5 13:39:04 2023                  *
***************************************************************************
                           Input parameters                                
***************************************************************************"""
    assert _read_sparc_version(header1) == "2023.02.03"
    header2 = """***************************************************************************
*                       SPARC (version June 24, 2023)                      *
*   Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech   *
*           Distributed under GNU General Public License 3 (GPL)          *
*                   Start time: Sun Feb  5 13:39:04 2023                  *
***************************************************************************
                           Input parameters                                
***************************************************************************"""
    assert _read_sparc_version(header2) == "2023.06.24"


def test_output_parser():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.out import _read_out

    data_dict = _read_out(
        test_output_dir
        / "AlSi_primitive_quick_relax.sparc"
        / "AlSi_primitive_quick_relax.out"
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


def test_output_parser_all():
    from sparc.common import repo_dir
    from sparc.sparc_parsers.out import _read_out

    for f_out in test_output_dir.glob("**/*.out"):
        data_dict = _read_out(f_out)
