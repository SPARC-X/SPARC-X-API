"""Test reading of all SPARC calculation examples

The test needs to be combined with the downloadable test outputs from SPARC-X's repo
and will only be activated when the environment variable $SPARC_TESTS_DIR is set

The ref
"""
import pytest
from pathlib import Path
import os
import tempfile
import shutil


def test_read_all_tests():
    """Search all .inpt files within the tests dir."""

    from sparc.io import read_sparc

    skipped_names = []
    tests_dir = os.environ.get("SPARC_TESTS_DIR", "")
    if len(tests_dir) == 0:
        pytest.skip(allow_module_level=True)

    tests_dir = Path(tests_dir)
    for inpt_file in tests_dir.glob("**/*.inpt"):
        workdir = inpt_file.parent
        parent_name = inpt_file.parents[1].name
        if parent_name in skipped_names:
            continue
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for ext in [".ion", ".inpt"]:
                shutil.copy(workdir / f"{parent_name}{ext}", tmpdir)
            for ext in [".refout", "refstatic", "refgeopt", "refaimd"]:
                origin_file = workdir / f"{parent_name}{ext}"
                new_ext = ext.replace("ref", "")
                new_file = tmpdir / f"{parent_name}{new_ext}"
                if origin_file.is_file():
                    shutil.copy(origin_file, new_file)
            try:
                read_sparc(tmpdir)
                # print("Passed: ", parent_name, workdir)
            except Exception as e:
                print("Failed: ", parent_name, workdir)
                print("\t: Error is ", e)
