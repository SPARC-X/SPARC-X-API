"""Convert ase atoms to structured dict following SPARC format
and vice versa
"""
import re

import numpy as np

from ase import Atoms, Atom
from ase.utils import reader, writer
from ase.io.utils import ImageIterator
from ase.io import ParseError
from ase.units import Bohr
from pathlib import Path

# from .sparc_parsers.ion import read_ion, write_ion

from ase import io
import numpy as np
from pathlib import Path
import shutil
import tarfile

from .ion import _read_ion, _write_ion, _ion_coord_to_ase_pos
from .inpt import _read_inpt, _write_inpt, _inpt_cell_to_ase_cell
from .utils import make_reverse_mapping
from ..inputs import SparcInputs

from warnings import warn



def atoms_to_dict(atoms, sort=True, direct=False, wrap=False, ignore_constraints=False, psp_dir=None, comments=""):
    """Given an ASE Atoms object, convert to SPARC ion and inpt data dict
    """
    # Step 1: if we should sort the atoms?
    origin_atoms = atoms.copy()
    if sort:
        sort_ = np.argsort(atoms.get_chemical_symbols())
        resort_ = make_reverse_mapping(sort_)
        # This is the sorted atoms object
        atoms = atoms[sort_]
    else:
        sort_ = []
        resort_ = []

    # Step 2: determine the counts of each element
    symbol_counts = count_symbols(atoms.get_chemical_symbols())
    write_spin = np.any(atoms.get_initial_magnetic_moments() != 0)
    has_charge = np.any(atoms.get_initial_charges() != 0)
    if has_charge:
        warn(("SPARC currently doesn't support changing total number of electrons! "
              "via nomimal charges. The initial charges in the structure will be ignored."
              ))

    atom_blocks = []
    # Step 3: write each block
    for symbol, start, end in symbol_counts:
        block_dict = {}
        block_dict["ATOM_TYPE"] = symbol
        block_dict["N_TYPE_ATOM"] = end - start
        # TODO: make pseudo finding work
        block_dict["PSEUDO_POT"] = f"{symbol}.psp8"
        # TODO: atomic mass?
        p_atoms = atoms[start: end]
        if direct:
            pos = p_atoms.get_scaled_positions(wrap=wrap)
            block_dict["COORD_FRAC"] = pos
        else:
            # TODO: should we use default converter?
            pos = p_atoms.get_positions(wrap=wrap) / Bohr
            block_dict["COORD"] = pos
        if write_spin:
            # TODO: should we process atoms with already calculated magmoms?
            block_dict["SPIN"] = p_atoms.get_initial_magnetic_moments()
        # TODO: get write_relax
        atom_blocks.append(block_dict)

    # Step 4: inpt part
    # TODO: what if atoms does not have cell?
    cell_au = atoms.cell / Bohr
    inpt_blocks = {"LATVEC": cell_au, "LATVEC_SCALE": [1., 1., 1.]}

    comments = comments.split("\n")
    ion_data = {"ion_atom_blocks": atom_blocks, "ion_comments": comments,
                "sorting": {"sort": sort_, "resort": resort_}}
    inpt_data = {"inpt_blocks": inpt_blocks, "inpt_comments": []}
    return {**ion_data, **inpt_data}


def dict_to_atoms(data_dict):
    """Given a SPARC struct dict, construct the ASE atoms object

    Note: this method supports only 1 Atoms at a time
    """
    ase_cell = _inpt_cell_to_ase_cell(data_dict["inpt_blocks"])
    new_atom_blocks = _ion_coord_to_ase_pos(data_dict["ion_atom_blocks"], ase_cell)
    # Now the real thing to construct an atom object
    atoms = Atoms()
    atoms.cell = ase_cell
    for block in new_atom_blocks:
        element = block["ATOM_TYPE"]
        positions = block["_ase_positions"]
        if positions.ndim == 1:
            positions = positions.reshape(1, -1)
            # Consider moving spins to another function
        spins = block.get("SPIN", None)
        if spins is None:
            spins = np.zeros_like(positions)
        for pos, spin in zip(positions, spins):
            # TODO: What about charge?
            atoms.append(Atom(symbol=element, position=pos, magmom=spin))

    # Now we do a sort on the atom indices. The atom positions read from
    # .ion correspond to the `sort` and we use `resort` to transform

    # TODO: should we store the sorting information in SparcBundle?
    if "sorting" in data_dict:
        resort = data_dict["sorting"]["resort"]
    if len(resort) != len(atoms):
        # TODO: new exception
        raise ValueError("Length of resort mapping is different from the number of atoms!")
    atoms = atoms[resort]
    # TODO: set pbc and relax
    atoms.pbc = True
    return atoms


def count_symbols(symbols):
    """Count the number of consecutive elements.
    Output tuple is: element, start, end 
    For example, "CHCHHO" --> [('C', 0, 1), ('H', 1, 2), ('C', 2, 3), ('H', 3, 5), ('O', 5, 6)]
    """
    counts = []
    current_count = 1
    current_symbol = symbols[0]
    for i, symbol in enumerate(symbols[1:], start=1):
        if symbol == current_symbol:
            current_count += 1
        else:
            counts.append((current_symbol, i - current_count, i))
            current_count = 1
            current_symbol = symbol
    i += 1
    counts.append((current_symbol, i - current_count, i))
    return counts
