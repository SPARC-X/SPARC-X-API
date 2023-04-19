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
from ase.constraints import FixAtoms, FixedLine, FixedPlane, FixScaled

from warnings import warn


def atoms_to_dict(
    atoms,
    sort=True,
    direct=False,
    wrap=False,
    ignore_constraints=False,
    psp_dir=None,
    comments="",
):
    """Given an ASE Atoms object, convert to SPARC ion and inpt data dict"""
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
        warn(
            (
                "SPARC currently doesn't support changing total number of electrons! "
                "via nomimal charges. The initial charges in the structure will be ignored."
            )
        )

    relax_mask = relax_from_all_constraints(atoms.constraints, len(atoms))
    write_relax = (len(relax_mask) > 0) and (not ignore_constraints)

    atom_blocks = []
    # Step 3: write each block
    for symbol, start, end in symbol_counts:
        block_dict = {}
        block_dict["ATOM_TYPE"] = symbol
        block_dict["N_TYPE_ATOM"] = end - start
        # TODO: make pseudo finding work
        block_dict["PSEUDO_POT"] = f"{symbol}.psp8"
        # TODO: atomic mass?
        p_atoms = atoms[start:end]
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
        if write_relax:
            relax_this_block = relax_mask[start:end]
            block_dict["RELAX"] = relax_this_block
        # TODO: get write_relax
        atom_blocks.append(block_dict)

    # Step 4: inpt part
    # TODO: what if atoms does not have cell?
    cell_au = atoms.cell / Bohr
    inpt_blocks = {"LATVEC": cell_au, "LATVEC_SCALE": [1.0, 1.0, 1.0]}

    comments = comments.split("\n")
    ion_data = {
        "ion_atom_blocks": atom_blocks,
        "ion_comments": comments,
        "sorting": {"sort": sort_, "resort": resort_},
    }
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
    relax_dict = {}

    atoms_count = 0
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
        relax = block.get("RELAX", [])
        for i, r in enumerate(relax, start=atoms_count):
            relax_dict[i] = r
        atoms_count += len(positions)

    if "sorting" in data_dict:
        resort = data_dict["sorting"]["resort"]
    else:
        resort = np.arange(len(atoms))

    if len(resort) != len(atoms):
        # TODO: new exception
        raise ValueError(
            "Length of resort mapping is different from the number of atoms!"
        )
    # TODO: check if this mapping is correct
    resorted_relax_dict = {resort[i]: r for i, r in relax_dict.items()}
    # Now we do a sort on the atom indices. The atom positions read from
    # .ion correspond to the `sort` and we use `resort` to transform

    # TODO: should we store the sorting information in SparcBundle?

    atoms = atoms[resort]
    constraints = constraints_from_relax(resorted_relax_dict)
    atoms.constraints = constraints

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


def constraints_from_relax(relax_dict):
    """
    Convert the SPARC RELAX fields to ASE's constraints

    Arguments
    relax: bool vector of size Nx3, i.e. [[True, True, True], [True, False, False]]

    Supported ase constraints will be FixAtoms, FixedLine and FixedPlane.
    For constraints in the same direction, all indices will be gathered.

    Note: ase>=3.22 will have FixedLine and FixedPlane accepting only 1 index at a time!

    The relax vector must be already sorted!
    """
    if len(relax_dict) == 0:
        return []

    cons_list = []
    # gathered_indices is an intermediate dict that contains
    # key: relax mask if not all True
    # indices: indices that share the same mask
    #
    gathered_indices = {}

    for i, r in relax_dict.items():
        r = tuple(np.ndarray.tolist(r.astype(bool)))
        if np.all(r):
            continue

        if r not in gathered_indices:
            gathered_indices[r] = [i]
        else:
            gathered_indices[r].append(i)

    for relax_type, indices in gathered_indices.items():
        degree_freedom = 3 - relax_type.count(True)

        if degree_freedom == 0:
            cons_list.append(FixAtoms(indices=indices))
        elif degree_freedom == 1:
            for ind in indices:
                cons_list.append(FixedLine(ind, np.array(relax_type).astype(int)))
        elif degree_freedom == 2:
            for ind in indices:
                cons_list.append(FixedPlane(ind, (~np.array(relax_type)).astype(int)))
    return cons_list


def relax_from_constraint(constraint):
    """returns dict of {atom_index: relax_dimensions} for the given constraint"""
    type_name = constraint.todict()["name"]
    if isinstance(constraint, FixAtoms):
        dimensions = [False] * 3
        expected_free = 0
    elif isinstance(constraint, FixedLine):
        # Only supports orthogonal basis!
        dimensions = [d == 1 for d in constraint.dir]
        expected_free = 1
    elif isinstance(constraint, FixedPlane):
        dimensions = [d != 1 for d in constraint.dir]
        expected_free = 2
    else:
        warn(
            f"The constraint type {type_name} is not supported by"
            " SPARC's .ion format. This constraint will be"
            " ignored"
        )
        return {}
    if dimensions.count(True) != expected_free:
        warn(
            "SPARC's .ion filetype can only support freezing entire "
            f"dimensions (x,y,z). The {type_name} constraint will be ignored"
        )
    return {i: dimensions for i in constraint.get_indices()}  # atom indices


def relax_from_all_constraints(constraints, natoms):
    """converts ASE atom constraints to SPARC relaxed dimensions for the atoms"""
    if len(constraints) == 0:
        return []

    relax = [
        [True, True, True],
    ] * natoms  # assume relaxed in all dimensions for all atoms
    for c in constraints:
        for atom_index, rdims in relax_from_constraint(c).items():
            # There might be multiple constraints applied on one index,
            # always make it more constrained
            relax[atom_index] = np.bitwise_and(relax[atom_index], rdims)
    return relax
