import pytest
import numpy as np
from pathlib import Path
import os
import tempfile

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

    monkeypatch.setattr(
        sparc_io_bundle, "is_psp_download_complete", _fake_psp_check
    )

    from sparc.io import SparcBundle
    from sparc.common import psp_dir as default_psp_dir

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

    sb = SparcBundle(
        directory=test_output_dir / "Cu_FCC.sparc", label="Something"
    )
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
    from sparc.io import SparcBundle
    from ase.units import Bohr, Angstrom

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

    sb = SparcBundle(
        directory=test_output_dir / "AlSi_primitive_quick_relax.sparc"
    )
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "AlSi"

    sb = SparcBundle(directory=test_output_dir / "Fe2_spin_scan_gamma.sparc")
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "Fe2"
    assert tuple(atoms.get_initial_magnetic_moments()) == (1.0, 1.0)

    sb = SparcBundle(
        directory=test_output_dir / "TiO2_orthogonal_quick_md.sparc"
    )
    atoms = sb._read_ion_and_inpt()
    assert atoms.get_chemical_formula() == "O4Ti2"


def test_write_ion_inpt(fs):
    """Same example as in test_parse_atoms but try writing inpt and atoms"""
    from sparc.io import SparcBundle
    from ase.units import Bohr, Angstrom
    from ase.build import bulk

    fs.create_dir("test.sparc")
    # data_dict = {
    #     "ion": {
    #         "atom_blocks": [
    #             {
    #                 "ATOM_TYPE": "Cu",
    #                 "ATOMIC_MASS": "63.546",
    #                 "PSEUDO_POT": "../../../psps/29_Cu_19_1.7_1.9_pbe_n_v1.0.psp8",
    #                 "N_TYPE_ATOM": 4,
    #                 "COORD_FRAC": np.array(
    #                     [
    #                         [0.0, 0.0, 0.0],
    #                         [0.0, 0.5, 0.5],
    #                         [0.5, 0.0, 0.5],
    #                         [0.5, 0.5, 0.0],
    #                     ]
    #                 ),
    #                 "RELAX": np.array(
    #                     [
    #                         [True, True, True],
    #                         [False, False, False],
    #                         [False, False, False],
    #                         [False, False, False],
    #                     ]
    #                 ),
    #             },
    #         ],
    #         "comments": [
    #             "=========================",
    #             "format of ion file",
    #             "=========================",
    #             "ATOM_TYPE: <atom type name>",
    #             "N_TYPE_ATOM: <num of atoms of this type>",
    #             "COORD:",
    #             "<xcoord> <ycoord> <zcoord>",
    #             "...",
    #             "RELAX:",
    #             "<xrelax> <yrelax> <zrelax>",
    #             "...",
    #             "atom type",
    #             "atomic mass (amu)",
    #             "pseudopotential file",
    #             "number of atoms of this type",
    #             "COORD:                      # Cartesian coordinates (au)",
    #             "fractional coordinates (in lattice vector basis)",
    #         ],
    #         "sorting": {"sort": [3, 2, 1, 0], "resort": [3, 2, 1, 0]},
    #     },
    #     "inpt": {
    #         "params": {
    #             "LATVEC": [
    #                 [5.5, 0, 0],
    #                 [0, 5.5, 0],
    #                 [0, 0, 5.5],
    #             ],
    #             "ELEC_TEMP": 300.0,
    #             "MESH_SPACING": 0.22,
    #         },
    #         "comments": [],
    #     },
    # }
    atoms = bulk("Cu") * [4, 4, 4]
    with pytest.raises(ValueError):
        sp = SparcBundle(directory="test.sparc", mode="r")
        sp._write_ion_and_inpt(atoms)

    sp = SparcBundle(directory="test.sparc", mode="w")
    sp._write_ion_and_inpt(atoms, direct=True, copy_psp=False)
    # Copy psp should have the psps available


def test_write_ion_inpt_real():
    """Same example as in test_parse_atoms but try writing inpt and atoms"""
    from sparc.io import SparcBundle
    from ase.units import Bohr, Angstrom
    from ase.build import bulk

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
