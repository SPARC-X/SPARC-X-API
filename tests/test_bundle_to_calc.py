import pytest
import numpy as np
from pathlib import Path
import os
import tempfile

curdir = Path(__file__).parent
test_psp_dir = curdir / "psps"
test_output_dir = curdir / "outputs"


def test_bundle_convert():
    from sparc.io import SparcBundle

    for bundle in test_output_dir.glob("*.sparc"):
        sp = SparcBundle(directory=bundle)
        sp.convert_to_ase()
    return


def test_bundle_read():
    from sparc.io import read_sparc

    for bundle in test_output_dir.glob("*.sparc"):
        images = read_sparc(bundle)
    return
