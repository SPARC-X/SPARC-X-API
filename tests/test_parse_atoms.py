from pathlib import Path

import numpy as np
import pytest

curdir = Path(__file__).parent
test_output_dir = curdir / "outputs"


def test_atoms2dict():
    from ase.build import molecule
    from ase.constraints import FixAtoms

    from sparc.sparc_parsers.atoms import atoms_to_dict

    mol = molecule("C2H6", pbc=True, cell=[10, 10, 10])

    adict = atoms_to_dict(mol, sort=True, direct=False, comments=["Ethane"])
    assert adict["ion"]["atom_blocks"][0]["ATOM_TYPE"] == "C"
    assert "dummy" in adict["ion"]["atom_blocks"][0]["PSEUDO_POT"]
    assert adict["ion"]["atom_blocks"][1]["ATOM_TYPE"] == "H"
    assert "dummy" in adict["ion"]["atom_blocks"][1]["PSEUDO_POT"]

    adict = atoms_to_dict(mol, sort=False, direct=True, comments=["Ethane"])
    assert adict["ion"]["atom_blocks"][0]["ATOM_TYPE"] == "C"
    assert "dummy" in adict["ion"]["atom_blocks"][0]["PSEUDO_POT"]
    assert adict["ion"]["atom_blocks"][1]["ATOM_TYPE"] == "H"
    assert "dummy" in adict["ion"]["atom_blocks"][1]["PSEUDO_POT"]

    # Charge atoms
    mol1 = mol.copy()
    mol1.set_initial_charges([0.1] * 8)
    with pytest.warns(UserWarning, match="initial charges"):
        # Charge cannot be written
        adict = atoms_to_dict(mol1, sort=False, direct=True, comments=["Ethane"])

    # Spin
    mol2 = mol.copy()
    mol2.set_initial_magnetic_moments([1.0] * 8)
    adict = atoms_to_dict(mol2, sort=True, direct=True, comments=["Ethane"])
    assert np.isclose(
        adict["ion"]["atom_blocks"][0]["SPIN"], np.array([1.0, 1.0])
    ).all()

    # Fix
    mol3 = mol.copy()
    mol3.constraints = [
        FixAtoms(
            [
                0,
                1,
            ]
        )
    ]
    adict = atoms_to_dict(mol3, sort=True, direct=True, comments=["Ethane"])
    print(np.array(adict["ion"]["atom_blocks"][0]["RELAX"]))
    print(np.array(adict["ion"]["atom_blocks"][1]["RELAX"]))
    assert (~np.array(adict["ion"]["atom_blocks"][0]["RELAX"])).all()
    assert (np.array(adict["ion"]["atom_blocks"][1]["RELAX"])).all()

    # Special case, only 1 atom
    mol = molecule("CH4", pbc=True, cell=[10, 10, 10])
    atoms_to_dict(mol, sort=False, direct=True, comments=["Methane"])


def test_dict2atoms():
    from ase.build import molecule
    from ase.constraints import FixAtoms
    from ase.units import Angstrom, Bohr

    from sparc.sparc_parsers.atoms import dict_to_atoms

    data_dict = {
        "ion": {
            "atom_blocks": [
                {
                    "ATOM_TYPE": "Cu",
                    "ATOMIC_MASS": "63.546",
                    "PSEUDO_POT": "../../../psps/29_Cu_19_1.7_1.9_pbe_n_v1.0.psp8",
                    "N_TYPE_ATOM": 4,
                    "COORD_FRAC": np.array(
                        [
                            [0.0, 0.0, 0.0],
                            [0.0, 0.5, 0.5],
                            [0.5, 0.0, 0.5],
                            [0.5, 0.5, 0.0],
                        ]
                    ),
                    "RELAX": np.array(
                        [
                            [True, True, True],
                            [False, False, False],
                            [False, False, False],
                            [False, False, False],
                        ]
                    ),
                },
            ],
            "comments": [
                "=========================",
                "format of ion file",
                "=========================",
                "ATOM_TYPE: <atom type name>",
                "N_TYPE_ATOM: <num of atoms of this type>",
                "COORD:",
                "<xcoord> <ycoord> <zcoord>",
                "...",
                "RELAX:",
                "<xrelax> <yrelax> <zrelax>",
                "...",
                "atom type",
                "atomic mass (amu)",
                "pseudopotential file",
                "number of atoms of this type",
                "COORD:                      # Cartesian coordinates (au)",
                "fractional coordinates (in lattice vector basis)",
            ],
            "sorting": {"sort": [3, 2, 1, 0], "resort": [3, 2, 1, 0]},
        },
        "inpt": {
            "params": {
                "LATVEC": [
                    [5.5, 0, 0],
                    [0, 5.5, 0],
                    [0, 0, 5.5],
                ],
                "ELEC_TEMP": 300.0,
                "GRID_SPACING": 0.22,
            },
            "comments": [],
        },
    }
    atoms = dict_to_atoms(data_dict)

    # constraints
    assert len(atoms.constraints) == 1
    assert atoms.constraints[0].todict()["name"] == "FixAtoms"
    assert set(atoms.constraints[0].get_indices()) == set([0, 1, 2])

    # Reconver the indices
    expected_pos = (
        np.array(
            [
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
                [0.0, 0.0, 0.0],
            ]
        )
        * 5.5
        * Bohr
    )
    assert np.isclose(atoms.positions, expected_pos).all()


def test_count_symbols():
    from sparc.sparc_parsers.atoms import count_symbols

    symbols = list("CHCHHO")
    expected_output = [
        ("C", 0, 1),
        ("H", 1, 2),
        ("C", 2, 3),
        ("H", 3, 5),
        ("O", 5, 6),
    ]
    assert count_symbols(symbols) == expected_output

    symbols = ["He"]
    expected_output = [("He", 0, 1)]
    assert count_symbols(symbols) == expected_output

    symbols = list("CHHHH")
    expected_output = [("C", 0, 1), ("H", 1, 5)]
    assert count_symbols(symbols) == expected_output


def test_constraint_from_relax():
    from sparc.sparc_parsers.atoms import constraints_from_relax

    relax_dict = {
        0: np.array([True, True, True]),
        1: np.array([False, False, False]),
        2: np.array([False, False, True]),  # Relax along z --> Fixed line
        3: np.array([False, True, True]),
    }  # Relax along yz --> Fixed plane yz
    cons = constraints_from_relax(relax_dict)
    assert len(cons) == 3
    print("CONSTRAINTS: ", cons)

    assert len(constraints_from_relax({})) == 0


def test_relax_from_constraint():
    from ase.constraints import FixAtoms, FixedLine, FixedPlane, Hookean

    from sparc.sparc_parsers.atoms import (
        relax_from_all_constraints,
        relax_from_constraint,
    )

    cons = [FixAtoms([0, 1, 2, 3])]
    assert (~np.array(relax_from_all_constraints(cons, 4))).all()
    with pytest.raises(ValueError):
        relax_from_all_constraints(cons, 3)

    assert tuple(relax_from_constraint(FixedLine(0, [0, 0, 1]))[0]) == (
        False,
        False,
        True,
    )
    assert tuple(relax_from_constraint(FixedPlane(0, [0, 0, 1]))[0]) == (
        True,
        True,
        False,
    )

    with pytest.warns(UserWarning, match="not supported"):
        assert relax_from_constraint(Hookean(0, 1, 10.0, 2.0)) == {}

    with pytest.warns(UserWarning, match="only support freezing entire"):
        # Freezing diagonal is not possible
        assert relax_from_constraint(FixedLine(0, [1, 1, 1])) == {}

    with pytest.warns(UserWarning, match="only support freezing entire"):
        # Freezing diagonal is not possible
        assert relax_from_constraint(FixedPlane(0, [1, 1, 1])) == {}


def test_atoms_pbc_conversion():
    from ase.build import bulk, molecule, mx2

    from sparc.sparc_parsers.atoms import atoms_to_dict

    h2 = molecule("H2")
    sparc_dict = atoms_to_dict(h2)

    assert sparc_dict["inpt"]["params"]["BC"] == "D D D"
    h2.pbc = True

    sparc_dict = atoms_to_dict(h2)
    assert sparc_dict["inpt"]["params"]["BC"] == "P P P"

    al = bulk("Al", cubic=True)
    sparc_dict = atoms_to_dict(al)
    assert sparc_dict["inpt"]["params"]["BC"] == "P P P"

    mos2 = mx2("MoS2")
    sparc_dict = atoms_to_dict(mos2)
    assert sparc_dict["inpt"]["params"]["BC"] == "P P D"


def test_atoms_sorting_order():
    """Chemical symbols should be sorted in a 'stable'-fashion
    i.e. the relative ordering between elements of the same symbol
    should appear the same after the sorting
    """
    from ase.build import bulk, molecule

    from sparc.sparc_parsers.atoms import atoms_to_dict

    # Case 1: All the same atoms, sort and resort should be the original order
    # stable-sort should retain the original order
    cu = bulk("Cu") * [5, 5, 5]
    ion_data = atoms_to_dict(cu, sort=True)["ion"]
    assert tuple(ion_data["sorting"]["sort"]) == tuple(range(125))
    assert tuple(ion_data["sorting"]["resort"]) == tuple(range(125))

    # Case 2: after sorting, the relative order of the same element should be incremental
    h2o = molecule("H2O", cell=[4, 4, 4], pbc=True) * [5, 5, 5]
    # The first 125 * 2 = 250 --> H
    ion_data = atoms_to_dict(h2o, sort=True)["ion"]
    sort_H = np.array(ion_data["sorting"]["sort"][:250])
    sort_O = np.array(ion_data["sorting"]["sort"][250:])
    # sort must be strictly incremental
    assert np.all(sort_H[1:] - sort_H[:-1])
    assert np.all(sort_O[1:] - sort_O[:-1])
