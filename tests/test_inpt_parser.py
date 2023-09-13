import pytest
from pathlib import Path
import numpy as np

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_read_inpt():
    from sparc.sparc_parsers.inpt import _read_inpt

    data_dict = _read_inpt(test_output_dir / "Cu_FCC.sparc/Cu_FCC.inpt")
    assert "inpt" in data_dict
    assert "params" in data_dict["inpt"]
    inpt_dict = data_dict["inpt"]
    params = inpt_dict["params"]
    assert np.isclose(params["LATVEC_SCALE"], np.ones(3) * 5.416914).all()
    assert np.isclose(params["LATVEC"], np.eye(3)).all()
    assert params["EXCHANGE_CORRELATION"] == "GGA_PBE"
    assert params["ELEC_TEMP_TYPE"] == "Fermi-Dirac"
    # TODO: BC may subject to change
    assert params["BC"] == "P P P"
    assert np.isclose(params["KPOINT_GRID"], np.ones(3)).all()

    assert params["MESH_SPACING"] == 0.22

    assert params["PRECOND_KERKER_THRESH"] == 0.0
    assert params["MIXING_PARAMETER"] == 1.0
    assert isinstance(params["MIXING_PARAMETER"], float)
    assert params["TOL_SCF"] == 0.0001
    assert params["CALC_STRESS"] is False
    assert params["PRINT_ATOMS"] is True
    assert params["PRINT_VELS"] is False
    assert params["PRINT_FORCES"] is False
    assert params["PRINT_EIGEN"] is False
    assert params["PRINT_DENSITY"] is False


def test_write_inpt():
    from sparc.sparc_parsers.inpt import _read_inpt, _write_inpt
    import tempfile

    data_dict = _read_inpt(test_output_dir / "Cu_FCC.sparc/Cu_FCC.inpt")
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        temp_inpt = tempdir / "temp.inpt"

        with pytest.raises(ValueError):
            _write_inpt(temp_inpt, {})

        with pytest.raises(ValueError):
            _write_inpt(temp_inpt, {"inpt": {}})

        _write_inpt(temp_inpt, data_dict)
        new_data_dict = _read_inpt(temp_inpt)
        for key in (
            "BC",
            "EXCHANGE_CORRELATION",
            "ELEC_TEMP_TYPE",
            "MESH_SPACING",
            "PRECOND_KERKER_THRESH",
            "MIXING_PARAMETER",
            "TOL_SCF",
            "CALC_STRESS",
            "PRINT_ATOMS",
            "PRINT_VELS",
            "PRINT_FORCES",
            "PRINT_EIGEN",
            "PRINT_DENSITY",
        ):
            assert (
                data_dict["inpt"]["params"][key]
                == new_data_dict["inpt"]["params"][key]
            )
        for key in ("LATVEC", "LATVEC_SCALE"):
            assert np.isclose(
                data_dict["inpt"]["params"][key],
                new_data_dict["inpt"]["params"][key],
            ).all()


def test_cell_conversion():
    from sparc.sparc_parsers.inpt import _inpt_cell_to_ase_cell
    from ase.units import Bohr, Angstrom

    # 0. invalid
    # 1. valid, equivalent to LATVEC = diag(1.5, 1.5, 1.5)
    data_dict = {
        "inpt": {
            "params": {
                "CELL": [1.5, 1.5, 1.5],
            }
        }
    }
    cell = _inpt_cell_to_ase_cell(data_dict)
    assert np.isclose(cell, np.eye(3) * 1.5 * Bohr).all()

    # 1. invalid
    data_dict = {
        "inpt": {
            "params": {"CELL": [1.5, 1.5, 1.5], "LATVEC_SCALE": [1.5, 1.5, 1.5]}
        }
    }
    with pytest.raises(ValueError):
        _inpt_cell_to_ase_cell(data_dict)

    # 2: Only LATVEC
    data_dict = {
        "inpt": {"params": {"LATVEC": [[1.5, 0, 0], [0, 1.5, 0], [0, 0, 1.5]]}}
    }
    cell = _inpt_cell_to_ase_cell(data_dict)
    assert np.isclose(cell, np.eye(3) * 1.5 * Bohr).all()

    # 3. LATVEC + LATVEC_SCALE
    data_dict = {
        "inpt": {
            "params": {
                "LATVEC": [[1.5, 0, 0], [0, 1.5, 0], [0, 0, 1.5]],
                "LATVEC_SCALE": [2, 3, 1],
            }
        }
    }
    cell = _inpt_cell_to_ase_cell(data_dict)
    assert np.isclose(cell, np.diag([2, 3, 1]) * 1.5 * Bohr).all()

    # 4. LATVEC + CELL
    data_dict = {
        "inpt": {
            "params": {
                "LATVEC": [[1.5, 0, 0], [0, 1.5, 0], [0, 0, 1.5]],
                "CELL": [8, 8, 9],
            }
        }
    }
    cell = _inpt_cell_to_ase_cell(data_dict)

    assert np.isclose(cell, np.diag([8, 8, 9]) * Bohr).all()
