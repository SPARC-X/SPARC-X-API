import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

curdir = Path(__file__).parent
test_psp_dir = curdir / "psps"
test_output_dir = curdir / "outputs"


def test_bundle_psp(monkeypatch):
    """Test PSP settings"""

    # Disable the default psp searching mechanism
    from sparc import io as sparc_io_bundle

    monkeypatch.setattr(sparc_io_bundle, "default_psp_dir", "/tmp")

    from sparc.io import SparcBundle

    os.environ.pop("SPARC_PP_PATH", None)
    os.environ.pop("SPARC_PSP_PATH", None)
    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")

    with pytest.warns(UserWarning, match="re-download"):
        sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc", mode="w")
    assert sb.psp_dir is None

    os.environ["SPARC_PP_PATH"] = test_psp_dir.as_posix()
    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")
    assert sb.psp_dir.resolve() == test_psp_dir.resolve()

    # SPARC_PSP_PATH has higher priority
    os.environ["SPARC_PSP_PATH"] = test_psp_dir.parent.as_posix()
    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")
    assert sb.psp_dir.resolve() == test_psp_dir.parent.resolve()

    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc", psp_dir="./")
    assert sb.psp_dir.resolve() == Path(".").resolve()


def test_default_psp(monkeypatch):
    """Test if default location of psp are correct"""
    from sparc import io as sparc_io_bundle

    # Make the psp downloader check always true
    def _fake_psp_check(directory):
        return True

    monkeypatch.setattr(sparc_io_bundle, "is_psp_download_complete", _fake_psp_check)

    from sparc.common import psp_dir as default_psp_dir
    from sparc.io import SparcBundle

    os.environ.pop("SPARC_PP_PATH", None)
    os.environ.pop("SPARC_PSP_PATH", None)
    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")
    assert Path(sb.psp_dir).resolve() == Path(default_psp_dir).resolve()


def test_bundle_label():
    """Test bundle label"""
    from sparc.io import SparcBundle

    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")
    sb.label == "Cu_FCC"
    assert sb._indir(".ion").name == "Cu_FCC.ion"

    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc", label="Something")
    sb.label == "Something"

    assert sb._indir(ext=".ion").name == "Something.ion"
    assert sb._indir(ext="ion").name == "Something.ion"

    with pytest.warns(UserWarning, match="illegal characters"):
        sb = SparcBundle(
            directory=test_output_dir / "Cu_FCC.sparc", label="Something?L"
        )
    assert sb.label == "SPARC"


def test_read_ion_inpt():
    """Test ion and inpt read"""
    from ase.units import Angstrom, Bohr

    from sparc.io import SparcBundle

    sb = SparcBundle(directory=test_output_dir / "Cu_FCC.sparc")
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "Cu4"
    assert atoms.cell.cellpar()[0] == 5.416914 * Bohr
    # Already resorted
    assert tuple(atoms.constraints[0].get_indices()) == (3,)
    assert np.isclose(
        atoms.positions,
        np.array(
            [
                [
                    0.5,
                    0.5,
                    0.0,
                ],
                [
                    0.5,
                    0.0,
                    0.5,
                ],
                [
                    0.0,
                    0.5,
                    0.5,
                ],
                [0.0, 0.0, 0.0],
            ]
        )
        * 5.416914
        * Bohr,
    ).all()

    sb = SparcBundle(directory=test_output_dir / "AlSi_primitive_quick_relax.sparc")
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "AlSi"

    sb = SparcBundle(directory=test_output_dir / "Fe2_spin_scan_gamma.sparc")
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "Fe2"
    assert tuple(atoms.get_initial_magnetic_moments()) == (1.0, 1.0)

    sb = SparcBundle(directory=test_output_dir / "TiO2_orthogonal_quick_md.sparc")
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "O4Ti2"


def test_write_ion_inpt(fs):
    """Same example as in test_parse_atoms but try writing inpt and atoms"""
    from ase.build import bulk
    from ase.units import Angstrom, Bohr

    from sparc.io import SparcBundle

    fs.create_dir("test.sparc")
    atoms = bulk("Cu") * [4, 4, 4]
    with pytest.raises(ValueError):
        sp = SparcBundle(directory="test.sparc", mode="r")
        sp._write_ion_and_inpt(atoms)

    sp = SparcBundle(directory="test.sparc", mode="w")
    sp._write_ion_and_inpt(atoms, direct=True, copy_psp=False)
    # Copy psp should have the psps available


def test_write_ion_inpt_real():
    """Same example as in test_parse_atoms but try writing inpt and atoms"""
    from ase.build import bulk
    from ase.units import Angstrom, Bohr

    from sparc.io import SparcBundle

    # Even without SPARC_PP_PATH, the psp files should exist
    os.environ.pop("SPARC_PP_PATH", None)
    os.environ.pop("SPARC_PSP_PATH", None)

    atoms = bulk("Cu") * [4, 4, 4]
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        workdir = tmpdir / "test.sparc"
        sp = SparcBundle(directory=workdir, mode="w")
        sp._write_ion_and_inpt(atoms, direct=True, copy_psp=True)
        # Copy psp should have the psps available
        assert len(list(Path(workdir).glob("*.psp8"))) == 1


def test_bundle_diff_label(fs):
    from ase.build import bulk
    from ase.units import Angstrom, Bohr

    from sparc.io import SparcBundle

    fs.create_dir("test.sparc")
    atoms = bulk("Cu") * [4, 4, 4]
    sp = SparcBundle(directory="test.sparc", mode="w", label="Cu")
    sp._write_ion_and_inpt(atoms)

    sp = SparcBundle(directory="test.sparc", mode="r")
    assert sp.label == "Cu"
    assert sp.directory.name == "test.sparc"

    # Write another into the bundle
    sp = SparcBundle(directory="test.sparc", mode="w", label="Cu1")
    sp._write_ion_and_inpt(atoms)

    with pytest.raises(Exception):
        sp = SparcBundle(directory="test.sparc", mode="r")


def test_bundle_write_multi(fs):
    import numpy as np
    from ase.build import bulk
    from ase.io import read, write

    from sparc.io import read_sparc, write_sparc

    fs.create_dir("test.sparc")
    atoms = bulk("Cu") * [4, 4, 4]
    images = [atoms]
    write_sparc("test.sparc", atoms)
    write_sparc("test.sparc", images)
    write("test.sparc", atoms, format="sparc")
    write("test.sparc", images, format="sparc")
    images.append(atoms.copy())
    with pytest.raises(Exception):
        write_sparc("test.sparc", images)

    with pytest.raises(Exception):
        write("test.sparc", images, format="sparc")

    atoms2 = read_sparc("test.sparc")
    atoms2 = read("test.sparc", format="sparc")
    assert np.isclose(atoms.positions, atoms2.positions).all()


def test_bundle_psp():
    from sparc.io import SparcBundle

    for f in test_output_dir.glob("*.sparc"):
        sb = SparcBundle(f)
        sb.read_raw_results()
        assert len(sb.psp_data) > 0
        # Embedded psp8 files should have the psp data available
        if len(list(f.glob("*.psp8"))) > 0:
            for elem, psp_data in sb.psp_data.items():
                assert "symbol" in psp_data


def test_bundle_reuse_sorting():
    """sort=True should reuse the sorting information from bundle"""
    import tempfile
    from pathlib import Path

    from sparc.io import SparcBundle

    sb = SparcBundle(test_output_dir / "Cu_FCC.sparc")
    init_atoms = sb.convert_to_ase()

    atoms2 = init_atoms.copy()
    atoms2.rattle()
    sort1 = sb.sorting["sort"]

    # Hard hack to write another directory, though not really done manually
    with tempfile.TemporaryDirectory() as tmpdir:
        sb.directory = Path(tmpdir)
        sb.mode = "w"
        sb._write_ion_and_inpt(atoms=atoms2)
        # Create another bundle reader
        sb2 = SparcBundle(tmpdir)
        atoms2 = sb2.convert_to_ase()
        assert np.isclose(sort1, sb2.sorting["sort"]).all()


def test_bundle_nh3():
    """Sorted & constrained NH3

    Order of atoms in SPARC ion file is H, H, H, N
    Make sure the constraint is aso correct
    """
    import tempfile
    from pathlib import Path

    from ase.build import molecule
    from ase.constraints import FixAtoms

    from sparc.io import SparcBundle

    nh3 = molecule("NH3", cell=(6, 6, 6))
    nh3.constraints = [FixAtoms([0])]
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        sb = SparcBundle(directory=tmpdir, mode="w")
        sb._write_ion_and_inpt(atoms=nh3)
        # Read back, the sorting should remain the same
        sb2 = SparcBundle(directory=tmpdir, mode="r")
        atoms = sb2.convert_to_ase()
        assert tuple(sb2.sorting["resort"]) == (3, 0, 1, 2)
        assert tuple(atoms.get_chemical_symbols()) == ("N", "H", "H", "H")
        assert atoms.constraints[0].index[0] == 0


def test_bundle_geopt_low_dim_stress():
    """Test if the low-dimensional stress from geopt and out
    are consistent
    """
    from sparc.io import SparcBundle

    sb = SparcBundle(test_output_dir / "Alloy_geopt_ppd_bc.sparc")
    raw_ionic_steps = sb.read_raw_results()["out"]["ionic_steps"]
    # Max equivalent stress from .out file
    max_stress_equiv_out = np.array(
        [
            entry["maximum stress equiv. to periodic"]["value"]
            for entry in raw_ionic_steps
        ]
    )
    # Max equivalent stress from the readout atoms
    images = sb.convert_to_ase(":")
    max_stress_equiv_geopt = np.array(
        [np.abs(atoms.calc.results["stress_equiv"]).max() for atoms in images]
    )
    print("Max stress (eV/Ang^3) from .out file:", max_stress_equiv_out)
    print("Max stress (eV/Ang^3) from .geopt file:", max_stress_equiv_geopt)
    assert np.isclose(max_stress_equiv_geopt, max_stress_equiv_geopt, 1e-6).all()
