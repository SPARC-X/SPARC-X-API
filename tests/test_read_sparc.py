import pytest
import numpy as np
from pathlib import Path
import os

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"

def test_read_sparc_all():
    from sparc.io import read_sparc
    from sparc.common import repo_dir

    for bundle in test_output_dir.glob("*.sparc/"):
        if bundle.name not in ("Al_multi_geopt.sparc",):
            results = read_sparc(bundle)
        else:
            results = read_sparc(bundle, include_all_files=True)