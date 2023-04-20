import pytest


def test_atoms2dict():
    from sparc.sparc_parsers.atoms import atoms_to_dict
    from ase.build import molecule

    mol = molecule("C2H6", pbc=True, cell=[10, 10, 10])

    # adict = atoms_to_dict(mol, sort=True, direct=False, comments=["Ethane"])
