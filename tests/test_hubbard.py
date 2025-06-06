import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

curdir = Path(__file__).parent
repo_dir = curdir.parent
test_output_dir = curdir / "outputs"


def test_parse_multiple_hubbard_atoms():
    """_read_ion should return a data dict
    containing multiple tuples of atomic U-correction
    values
    """
    from sparc.sparc_parsers.ion import _read_ion, _write_ion

    ion_content = """
#=========================
# format of ion file
#=========================
# ATOM_TYPE: <atom type name>
# PSEUDO_POT: <path/to/pseudopotential/file>
# N_TYPE_ATOM: <num of atoms of this type>
# COORD:
# <xcoord> <ycoord> <zcoord>
# ...
# RELAX:
# <xrelax> <yrelax> <zrelax>
# ...

# Reminder: when changing number of atoms, change the RELAX flags accordingly
#           as well.

ATOM_TYPE: Ti                # atom type
PSEUDO_POT: ../../../psps/22_Ti_12_2.0_2.8_pbe_n_v1.0.psp8
N_TYPE_ATOM: 2               # number of atoms of this type
COORD_FRAC:                  # coordinates follows
0.000000 0.500000 0.250000
0.500000 0.000000 0.750000

ATOM_TYPE: Cr                # atom type
PSEUDO_POT: ../../../psps/24_Cr_14_1.7_2.1_pbe_n_v1.0.psp8
N_TYPE_ATOM: 2               # number of atoms of this type
COORD_FRAC:                  # coordinates follows
0.000000 0.000000 0.000000
0.500000 0.500000 0.500000

ATOM_TYPE: O               # atom type
PSEUDO_POT: ../../../psps/08_O_6_1.2_1.4_pbe_n_v1.0.psp8
N_TYPE_ATOM: 8               # number of atoms of this type
COORD_FRAC:                  # coordinates follows
0.000000 0.500000 0.460598
0.500000 0.500000 0.297638
0.000000 0.500000 0.039402
0.000000 0.000000 0.202362
0.500000 0.000000 0.960598
0.000000 0.000000 0.797638
0.500000 0.000000 0.539402
0.500000 0.500000 0.702362

HUBBARD:
#U_ATOM_TYPE: Ti
#U_VAL: 0 0 0.1 0 # U values in hartrees for s p d f orbitals (valence only)

U_ATOM_TYPE: Cr
U_VAL: 0 0 0.2 0 # U values in hartrees for s p d f orbitals (valence only)

U_ATOM_TYPE: Ti
U_VAL: 0 0 0.1 0 # U values in hartrees for s p d f orbitals (valence only)
    """
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        test_ion = tempdir / "test.ion"
        test_ion_w = tempdir / "test_w.ion"
        with open(test_ion, "w") as fd:
            fd.write(ion_content)

        data_dict = _read_ion(test_ion)
        ion_dict = data_dict["ion"]
        assert "extra" in ion_dict
        assert "hubbard" in ion_dict["extra"]
        hubbard_u_pairs = ion_dict["extra"]["hubbard"]
        assert len(hubbard_u_pairs) == 2
        # Atom types are read in sequence
        # comments should be ommitted
        assert hubbard_u_pairs[0]["U_ATOM_TYPE"] == "Cr"
        assert hubbard_u_pairs[1]["U_ATOM_TYPE"] == "Ti"
        # HUBBARD U_VAL should be 4-array (s,p,d,f)
        for pair in hubbard_u_pairs:
            assert len(pair["U_VAL"]) == 4

        # Write the ion file from raw data_dict
        _write_ion(test_ion_w, data_dict)
        with open(test_ion_w, "r") as fd:
            new_ion_content = [line for line in fd.readlines()]
        assert any(["HUBBARD:" in line for line in new_ion_content])

        u_atom_type_found = False
        for line in new_ion_content:
            if "U_ATOM_TYPE:" in line:
                u_atom_type_found = True
            if "U_VAL:" in line:
                if u_atom_type_found is False:
                    raise ValueError("U_VAL line comes before U_ATOM_TYPE!")
                u_atom_type_found = False


def test_ion_to_atoms():
    """Atoms objects returned from parser
    should contain info section about the HUBBARD parameters
    """
    from sparc.io import read_sparc

    atoms = read_sparc(test_output_dir / "MoO3_hubbard.sparc")
    print(atoms.info)
    assert "hubbard_u (hartree)" in atoms.info
    hubbard_info = atoms.info["hubbard_u (hartree)"]
    assert isinstance(hubbard_info, list)
    assert len(hubbard_info) == 1
    assert hubbard_info[0]["U_ATOM_TYPE"] == "Mo"
    assert len(hubbard_info[0]["U_VAL"]) == 4
    # The U-values from atoms.info should be using Hartree unit
    assert np.isclose(hubbard_info[0]["U_VAL"][2], 0.05)


def test_write_atoms_to_ion():
    """Write atoms objects that carries HUBBARD info
    to ion.
    Without setting HUBBARD_FLAG the HUBBARD blocks won't be
    written
    """
    from sparc.io import read_sparc, write_sparc

    atoms = read_sparc(test_output_dir / "MoO3_hubbard.sparc")
    print(atoms.info)
    assert "hubbard_u (hartree)" in atoms.info
    # write without HUBBARD
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        write_sparc(tmpdir, atoms, label="SPARC")
        with open(tmpdir / "SPARC.ion", "r") as fd:
            ion_content = fd.readlines()
        print(ion_content)
        assert all(["HUBBARD:" not in line for line in ion_content])

    # write with HUBBARD
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        write_sparc(tmpdir, atoms, label="SPARC", input_parameters={"HUBBARD_FLAG": 1})
        with open(tmpdir / "SPARC.ion", "r") as fd:
            ion_content = fd.readlines()
        print(ion_content)
        with open(tmpdir / "SPARC.inpt", "r") as fd:
            inpt_content = fd.readlines()
        print(inpt_content)
        assert ["HUBBARD:" in line for line in ion_content]
        assert ["U_ATOM_TYPE: Mo" in line for line in ion_content]


def test_calc_hubbard_block():
    """Allow write calculation using hubbard block"""
    from ase.units import Hartree

    from sparc.calculator import SPARC
    from sparc.io import read_sparc

    # We will create a input file with user-provided HUBBARD
    # blocks
    atoms = read_sparc(test_output_dir / "MoO3_hubbard.sparc")
    print(atoms.info["hubbard_u (hartree)"][0])
    input_params = {
        "HUBBARD": [{"U_ATOM_TYPE": "Mo", "U_VAL": [0, 0, 1.0, 0]}],
        "HUBBARD_FLAG": 1,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(directory=tmpdir, **input_params)
        calc.write_input(atoms)
        with open(tmpdir / "SPARC.ion", "r") as fd:
            ion_content = fd.readlines()
        assert any(["HUBBARD" in line for line in ion_content])
        atoms_read = read_sparc(tmpdir)
        assert "hubbard_u (hartree)" in atoms_read.info
        assert np.isclose(atoms_read.info["hubbard_u (hartree)"][0]["U_VAL"][2], 1.0)


def test_bad_hubbard_block_write():
    """Allow write calculation using hubbard block"""
    from ase.units import Hartree

    from sparc.calculator import SPARC
    from sparc.io import read_sparc

    # We will create a input file with user-provided HUBBARD
    # blocks
    atoms = read_sparc(test_output_dir / "MoO3_hubbard.sparc")
    atoms.info.pop("hubbard_u (hartree)", None)
    # Generate an error when writing .ion file for bad element
    input_params = {
        "HUBBARD": [{"U_ATOM_TYPE": "Ni", "U_VAL": [0, 0, 1.0, 0]}],
        "HUBBARD_FLAG": 1,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(directory=tmpdir, **input_params)
        with pytest.raises(ValueError):
            calc.write_input(atoms)

    # Generate an error when writing .ion file for duplicated element
    input_params = {
        "HUBBARD": [
            {"U_ATOM_TYPE": "Ni", "U_VAL": [0, 0, 1.0, 0]},
            {"U_ATOM_TYPE": "Ni", "U_VAL": [0, 0, 2.0, 0]},
        ],
        "HUBBARD_FLAG": 1,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(directory=tmpdir, **input_params)
        with pytest.raises(ValueError):
            calc.write_input(atoms)

    # Generate an error when writing .ion file for no hubbard info
    input_params = {
        "HUBBARD_FLAG": 1,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(directory=tmpdir, **input_params)
        with pytest.raises(ValueError):
            calc.write_input(atoms)


def test_gpaw_style_setups():
    """Check whether gpaw-style setups keyword can be written"""
    from ase.units import Hartree

    from sparc.calculator import SPARC
    from sparc.io import read_sparc

    atoms = read_sparc(test_output_dir / "MoO3_hubbard.sparc")
    # 5.0 eV HUBBARD U. For safety concern we will only write
    # HUBBARD settings when HUBBARD_FLAG=1
    # note the new HUBBARD U_VAL is overwritten!
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        calc = SPARC(setups={"Mo": ":d,5.0"}, hubbard_flag=1, directory=tmpdir)
        calc.write_input(atoms)
        with open(tmpdir / "SPARC.ion", "r") as fd:
            ion_content = fd.readlines()
        assert any(["HUBBARD" in line for line in ion_content])
        atoms_read = read_sparc(tmpdir)
        assert "hubbard_u (hartree)" in atoms_read.info
        print(atoms.info["hubbard_u (hartree)"][0])
        print(atoms_read.info["hubbard_u (hartree)"][0])
        # The u_val for Mo is 5.0 eV, now converted to Ha
        assert np.isclose(
            atoms_read.info["hubbard_u (hartree)"][0]["U_VAL"][2], 5.0 / Hartree
        )


def test_read_hubbard_raw_results():
    """Parse the outputfiles from HUBBARD calculations"""
    from ase.units import Hartree

    from sparc.io import SparcBundle, read_sparc

    sb = SparcBundle(test_output_dir / "MoO3_hubbard.sparc")
    raw_results = sb.read_raw_results()
    first_ionic_step = raw_results["out"]["ionic_steps"][0]
    assert "u correction" in first_ionic_step
    assert first_ionic_step["u correction"]["unit"] == "eV"
    assert np.isclose(
        first_ionic_step["u correction"]["value"], 8.6236340534e-02 * Hartree
    )
