"""Test reading of all SPARC calculation examples

The test needs to be combined with the downloadable test outputs from SPARC-X's repo
and will only be activated when the environment variable $SPARC_TESTS_DIR is set

The ref
"""
import pytest
import numpy as np
from pathlib import Path
import os
import tempfile
import shutil

skipped_names = ["Si2_domain_paral", "Si2_kpt_paral",
                 "SiH4", "SiH4_quick",
                 "H2O_sheet_quick", "H2O_sheet",
                 "CdS_bandstruct"]

def test_read_all_tests():
    """Search all .inpt files within the tests dir."""

    from sparc.io import read_sparc

    skipped_names = []
    tests_dir = os.environ.get("SPARC_TESTS_DIR", "")
    print(f"Current test dir is {tests_dir}")
    if len(tests_dir) == 0:
        pytest.skip(allow_module_level=True)

    tests_dir = Path(tests_dir)
    failed_counts = 0
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
                failed_counts += 1
    if failed_counts > 0:
        raise RuntimeError("More than 1 test in output read test failed")

def test_write_all_inputs():
    """Search all .inpt files within the tests dir."""

    from sparc.io import read_sparc, read_ion, write_ion
    from sparc.sparc_parsers.inpt import _read_inpt

    # Skipped tests are to avoid unwanted keywords
    tests_dir = os.environ.get("SPARC_TESTS_DIR", "")
    failed_counts = 0
    print(f"Current test dir is {tests_dir}")
    if len(tests_dir) == 0:
        pytest.skip(allow_module_level=True)

    tests_dir = Path(tests_dir)
    for inpt_file in tests_dir.glob("**/*.inpt"):
        workdir = inpt_file.parent
        parent_name = inpt_file.parents[1].name
        ion_file = inpt_file.with_suffix(".ion")
        if parent_name in skipped_names:
            continue
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            origin_atoms = read_ion(ion_file)
            origin_inpt_dict = _read_inpt(inpt_file)
            for key in ["CELL", "LATVEC_SCALE", "LATVEC", "BC"]:
                origin_inpt_dict["inpt"]["params"].pop(key, None)
            # Re-write the ion and inpt files
            try:
                write_ion(tmpdir / "test.ion", origin_atoms, **origin_inpt_dict["inpt"]["params"])
                new_atoms = read_ion(tmpdir / "test.ion")
                new_inpt_dict = _read_inpt(tmpdir / "test.inpt")
                assert np.all(origin_atoms.pbc == new_atoms.pbc)
                for key in origin_inpt_dict["inpt"]["params"].keys():
                    origin_val = origin_inpt_dict["inpt"]["params"][key]
                    new_val = new_inpt_dict["inpt"]["params"][key]
                    if isinstance(origin_val, (int, bool)):
                        assert origin_val == new_val
                    elif isinstance(origin_val, float):
                        assert np.isclose(origin_val, new_val, 1e-6)
                    elif isinstance(origin_val, str):
                        assert origin_val == new_val
                    # Vector types can be list compared
                    elif isinstance(origin_val, (list, np.ndarray)):
                        assert np.all(origin_val == new_val)
                    
            except Exception as e:
                print("Failed: ", parent_name, workdir)
                print("\t: Error is ", e)
                failed_counts += 1
    if failed_counts > 0:
        raise RuntimeError("More than 1 test in inpt write test failed")
