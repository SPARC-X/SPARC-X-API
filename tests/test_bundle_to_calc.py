import pytest
import numpy as np
from pathlib import Path
import os
import tempfile
from ase.units import Hartree

curdir = Path(__file__).parent
test_psp_dir = curdir / "psps"
test_output_dir = curdir / "outputs"


def test_bundle_convert():
    from sparc.io import SparcBundle

    for bundle in test_output_dir.glob("*.sparc"):
        if bundle.name not in ("Al_multi_geopt.sparc",):
            sp = SparcBundle(directory=bundle)
            sp.convert_to_ase()
        else:
            sp = SparcBundle(directory=bundle)
            sp.convert_to_ase(include_all_files=True)

    return


def test_bundle_read():
    from sparc.io import read_sparc

    for bundle in test_output_dir.glob("*.sparc"):
        if bundle.name not in ("Al_multi_geopt.sparc",):
            print(bundle)
            images = read_sparc(bundle)
        else:
            images = read_sparc(bundle, include_all_files=True)
    return


def test_multi_file_geopt_read():
    from sparc.io import read_sparc

    bundle = test_output_dir / "Al_multi_geopt.sparc"

    # Last image is empty, must use include_all_files=True
    with pytest.raises(Exception):
        images = read_sparc(bundle, include_all_files=False)

    # default
    with pytest.raises(Exception):
        images = read_sparc(bundle, include_all_files)

    images = read_sparc(bundle, index=":", include_all_files=True)
    assert len(images) == 7

    last = read_sparc(bundle, index=-1, include_all_files=True)
    assert np.isclose(last.get_potential_energy(), -9.057488887961474 * Hartree, 1e-4)
    return
