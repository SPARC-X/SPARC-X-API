import os
from pathlib import Path

import numpy as np
import pytest

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_read_sparc_all():
    from sparc.common import repo_dir
    from sparc.io import read_sparc

    for bundle in test_output_dir.glob("*.sparc/"):
        if bundle.name not in ("Al_multi_geopt.sparc",):
            results = read_sparc(bundle)
        else:
            results = read_sparc(bundle, include_all_files=True)


def test_atoms_read_pbc():
    from sparc.io import read_sparc
    from sparc.sparc_parsers.atoms import atoms_to_dict

    # Case 1: H2O sheet
    water_sheet = read_sparc(test_output_dir / "H2O_sheet_yz.sparc")
    assert all(water_sheet.pbc == [False, True, True])

    # Case 2: H2O wire
    water_wire = read_sparc(test_output_dir / "H2O_wire_z.sparc")
    assert all(water_wire.pbc == [False, False, True])
