import pytest
import numpy as np
from pathlib import Path
import os
from ase.units import Bohr, Hartree

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_static_parser(fs):
    from sparc.common import repo_dir

    # Pyfakefs requires following line to add the external data folders
    fs.add_real_directory(repo_dir)
    from sparc.sparc_parsers.static import _read_static

    fs.create_file("test.static")
    with open("test.static", "w") as fd:
        fd.write(
            """***************************************************************************
                            Atom positions                                 
***************************************************************************
Fractional coordinates of Fe:
      0.0000000000       0.0000000000       0.0000000000
      0.5200000000       0.5100000000       0.4900000000
Total free energy (Ha): -2.283157353113279E+02
Atomic forces (Ha/Bohr):
  8.0738249305E-01   3.7399117306E-01  -3.5796157735E-01
 -8.0738249305E-01  -3.7399117306E-01   3.5796157735E-01
Stress (GPa): 
 -2.1918863425E+04   1.3932450782E+03  -5.1023512490E+01 
  1.3932450782E+03  -2.1975897437E+04  -1.2676410947E+02 
 -5.1023512490E+01  -1.2676410947E+02  -2.2380745784E+04
        """
        )
    data_dict = _read_static("test.static")
    assert "static" in data_dict
    static_dict = data_dict["static"]
    assert "atoms" in static_dict
    assert tuple(static_dict["atoms"]["symbols"]) == ("Fe", "Fe")
    assert np.isclose(static_dict["free energy"], -2.283157353113279e02 * Hartree)
    assert static_dict["forces"].shape == (2, 3)
    assert np.isclose(static_dict["forces"][0, 0], 8.0738249305e-01 * Hartree / Bohr)
    assert static_dict["stress"].shape == (6,)


def test_static_parser_missing_fields(fs):
    from sparc.common import repo_dir

    # Pyfakefs requires following line to add the external data folders
    fs.add_real_directory(repo_dir)
    from sparc.sparc_parsers.static import _read_static, _add_cell_info

    fs.create_file("test.static")
    with open("test.static", "w") as fd:
        fd.write(
            """***************************************************************************
                            Atom positions                                 
***************************************************************************
Fractional coordinates of Al:
      0.0000000000       0.0000000000       0.0000000000
Total free energy (Ha): -2.212996080865029E+00
Atomic forces (Ha/Bohr):
  0.0000000000E+00   0.0000000000E+00   0.0000000000E+00
        """
        )
    data_dict = _read_static("test.static")
    assert "static" in data_dict
    static_dict = data_dict["static"]
    assert "atoms" in static_dict
    assert "stress" not in static_dict
    assert tuple(static_dict["atoms"]["symbols"]) == ("Al",)
    cell = (
        np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
        * 5.656854249492380
        * Bohr
    )
    new_data_dict = _add_cell_info(data_dict, cell)

    assert "coord" in new_data_dict["static"]["atoms"]
    assert new_data_dict["static"]["atoms"]["coord"].shape == (1, 3)
    assert np.isclose(new_data_dict["static"]["atoms"]["coord"][0, 0], 0)


def test_static_parser_no_atoms(fs):
    from sparc.common import repo_dir

    # Pyfakefs requires following line to add the external data folders
    fs.add_real_directory(repo_dir)
    from sparc.sparc_parsers.static import _read_static, _add_cell_info

    fs.create_file("test.static")
    with open("test.static", "w") as fd:
        fd.write(
            """
Total free energy (Ha): -2.212996080865029E+00
        """
        )
    data_dict = _read_static("test.static")
    assert "static" in data_dict
    static_dict = data_dict["static"]
    assert "atoms" not in static_dict
    assert "stress" not in static_dict
    cell = (
        np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
        * 5.656854249492380
        * Bohr
    )
    new_data_dict = _add_cell_info(data_dict, cell)

    assert "atoms" not in new_data_dict["static"]
