"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)

This file has been heavily modified since SPARC 0.1

TODO: more descriptions about this file io parser
"""
import shutil
import os
from typing import List
from collections import namedtuple
import warnings
from warnings import warn

import numpy as np
from ase import Atoms, Atom
from ase.units import Bohr
from ase.constraints import FixAtoms, FixedLine, FixedPlane

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from .utils import get_label, strip_comments, bisect_and_strip

from ..inputs import SparcInputs

defaultAPI = SparcInputs()

@reader
def _read_ion(fileobj):
    """
    Read information from the .ion file. Note, this method does not return an atoms object,
    but rather return a dict. Thus the label option is not necessary to keep


    Reads an ion file. Because some of the information necessary to create
    an atoms object is found in the .inpt file, this function also attemtps to read
    that as a source of data. If the file is not found or the information is invalid,
    it will look for it in the comments of the ion file, as written.
    """
    contents = fileobj.read()
    # label = get_label(fileobj, ".ion")
    data, comments = strip_comments(contents)
    # We do not read the cell at this time!
    
    # find the index for all atom type lines. They should be at the top of their block
    atom_type_bounds = [i for i, x in enumerate(data) if "ATOM_TYPE" in x] + [len(data)]
    atom_blocks = [
        read_atom_block(data[start:end])
        for start, end in zip(atom_type_bounds[:-1], atom_type_bounds[1:])
    ]
    
    return {"data": atom_blocks, "comments": comments}
    
    # try:
    #     comment_data = read_comments(comments)
    # except Exception as e:
    #     warnings.warn(f"failed to read comment data: {e}")
    #     comment_data = CommentData(cell=[], indices=[], pbc_list=[])
    # print(comment_data)

    # try:
    #     cell = read_inpt_cell(label + ".inpt")
    # except Exception as e:
    #     warnings.warn(f"Failed to read inpt file: {e}")
    #     cell = comment_data.cell

    # if len(cell) == 0:
    #     cell = np.zeros((3, 3))
    #     warnings.warn(
    #         "No lattice vectors were found in either the .inpt"
    #         " file or in the comments of the .ion file. Thus no"
    #         " unit cell was set for the resulting atoms object. "
    #         "Use the output atoms at your own risk"
    #     )

    # # add end of list to bound last block
    # atom_type_bounds.append(len(stripped))

    # # iterates over tuples (i_0, i_1), (i_1, i_2), ...

    # raw_atoms, relax, spins = process_atom_blocks(atom_blocks, cell)
    # print(raw_atoms, relax, spins)

    # check if we can reorganize the indices
    # if recover_indices and len(comment_data.indices) == len(raw_atoms):
    #     raw_atoms = reorder(raw_atoms, comment_data.indices)
    #     relax = reorder(relax, comment_data.indices)
    #     spins = reorder(spins, comment_data.indices)

    # atoms = Atoms()
    # atoms.cell = cell
    # if len(comment_data.pbc_list):
    #     atoms.set_pbc(comment_data.pbc_list)

    # for atom in raw_atoms:
    #     atoms.append(atom)

    # atoms.set_initial_magnetic_moments(spins)
    # if recover_constraints:
    #     constraints = constraints_from_relax(relax)
    #     atoms.set_constraint(constraints)
    # return atoms

@writer
def _write_ion(
        fileobj,
        data_dict,
):
    """
    Writes the ion file content from the atom_dict

    Please note this is not a Atoms-compatible function!

    The data_dict takes similar format as _read_ion

    Basically, we want to ensure
    data_dict = _read_ion("some.ion")
    _write_ion("some.ion", data_dict)
    shows the same format
    """
    if "data" not in data_dict:
        raise ValueError("Must provide a data-section in the data_dict (blocks of atomic information)")

    comments = data_dict.get("comments", "")
    if len(comments) == 0:
        comments = ["Input File Generated By SPARC ASE Calculator"]
    for line in comments:
        fileobj.write(f"# {line}\n")
    fileobj.write("\n")
    blocks = data_dict["data"]
    for block in blocks:
        for key in ["ATOM_TYPE", "N_TYPE_ATOM", "PSEUDO_POT", "COORD_FRAC", "COORD", "RELAX"]:
            val = block.get(key, None)
            print(key, val)
            if (key not in ["RELAX", "COORD", "COORD_FRAC"]) and (val is None):
                raise ValueError(f"Key {key} is not provided! Abort writing ion file")
            # TODO: change the API version
            if val is None:
                continue
            
            val_string = defaultAPI.convert_value_to_string(key, val)
            print(val_string)
            # TODO: make sure 1 line is accepted 
            if (val_string.count("\n") > 0) or (key in ["COORD_FRAC", "COORD", "RELAX"]):
                output = f"{key}:\n{val_string}\n"
            else:
                output = f"{key}: {val_string}\n"
            fileobj.write(output)
        # Write a split line
        fileobj.write("\n")
    return
    

    
    # env_pseudo = os.environ.get("SPARC_PSP_PATH", "")
    # for val in [pseudo_dir, env_pseudo]:
    #     if val and os.path.isdir(val):
    #         pseudo_dir = val
    #         break
    # else:
    #     # loop-else runs if loop doesn't break
    #     pseudo_dir = ""
    #     warnings.warn(
    #         "no valid value given for pseudo_dir. Explicitly set the argument or define the environment variable $SPARC_PSP_PATH to ensure your output references valid pseudopotential files"
    #     )

    # relax = (
    #     relax_from_all_constraints(atoms.constraints, len(atoms))
    #     if add_constraints
    #     else []
    # )

    # grouped_atoms = {}  # group by atom type
    # for atom in atoms:
    #     grouped_atoms.setdefault(atom.symbol, []).append(atom)

    # # header
    # fileobj.write("# Input File Generated By SPARC ASE Calculator #\n")
    # fileobj.write(format_cell(atoms.cell) + "\n")
    # fileobj.write(format_latvec(atoms.cell) + "\n")
    # if atoms.pbc is not None:
    #     fileobj.write(format_pbc(atoms.pbc) + " \n")
    # fileobj.write("# " + comment + "\n\n\n")

    # # body
    # directory = os.path.dirname(fileobj.name)
    # for element, block_atoms in sorted(grouped_atoms.items()):
    #     try:
    #         pseudo_path = find_pseudo_path(element, pseudo_dir)
    #         if copy_psp:
    #             copy_psp(pseudo_path, directory)
    #             pseudo_path = os.path.basename(pseudo_path)
    #     except:
    #         pseudo_path = element + ".pot"
    #         warnings.warn(
    #             "No pseudopotential detected for " + element + ". "
    #             "Using generic filename (" + pseudo_path + "). The pseudopotential "
    #             " must be added manually for SPARC to execute successfully."
    #         )
    #     block_string = format_atom_block(
    #         block_atoms, relax=relax, scaled=scaled, pseudo_path=pseudo_path
    #     )
    #     fileobj.write(f"{block_string}\n\n\n")


CommentData = namedtuple("CommentData", "cell indices pbc_list")


def read_comments(comments) -> CommentData:
    cell = []
    if "CELL" in comments[1]:  # Check only the second line for a CELL
        try:
            cell = comments[1].split()[1:]
            cell = [float(a) for a in cell]
            lat_array = read_lat_array(comments[3:])
            cell = (lat_array.T * cell).T * Bohr
        except:
            pass
    indices = []
    pbc_list = []
    for comment in comments:
        first, remaining = bisect_and_strip(comment, " ")
        if "index" == first:
            try:
                indices.append(int(remaining))
            except:
                pass
        if "PBC:" in comment:
            pbc_list = [c.strip().lower() == "true" for c in remaining.split()]

    return CommentData(cell=cell, indices=indices, pbc_list=pbc_list)


def read_atom_block(block, validator=defaultAPI):
    block_dict = {}
    multiline_key = ""
    use_validator = True if validator else False
    for line in block:
        if ":" not in line:
            # no key, assume multiline value
            block_dict[multiline_key].append(line.strip())
            continue
        key, value = bisect_and_strip(line, ":")
        key = key.upper()
        # print(key, value)
        if key and value:
            block_dict[key] = value
        elif key:
            # no value, assume that this key has a list of values
            # in the following lines
            block_dict[key] = []
            multiline_key = key
    validate_atom_block(block_dict)
    for key, val in block_dict.items():
        # print(key, val)
        _use_validator_this_key = use_validator
        if _use_validator_this_key:
            if key not in validator.parameters.keys():
                warn(f"Key {key} not in validator's parameter list, ignore value conversion!")
                _use_validator_this_key = False
        if _use_validator_this_key:
            val = validator.convert_string_to_value(key, val)
        block_dict[key] = val
    return block_dict


NativeData = namedtuple("NativeData", "atoms relax spins")


def process_atom_blocks(blocks, cell) -> NativeData:
    atoms = []
    relax = []
    spins = []
    for block in blocks:
        natoms = int(block["N_TYPE_ATOM"])
        if "COORD_FRAC" in block:
            for coords in block["COORD_FRAC"]:
                coords = coords.split()[:3]
                assert len(coords) == 3
                coords = sum([float(x) * a for x, a in zip(coords, cell)])
                atoms.append(Atom(symbol=block["ATOM_TYPE"], position=coords))
        elif "COORD" in block:
            for coords in block["COORD"]:
                coords = [float(x) * Bohr for x in block["COORD"].split()[:3]]
                assert len(coords) == 3
                atoms.append(Atom(symbol=block["ATOM_TYPE"], position=coords))
        if "RELAX" in block:
            for cons_set in block["RELAX"]:
                relax.append([int(a) for a in cons_set.split()])
        else:
            # there aren't constraints with this block, put in empty lists
            relax += [[]] * natoms
        if "SPIN" in block:
            spins += [float(spin) for spin in block["SPIN"]]
        else:
            # no spins in this block, set to 0
            spins += [0] * natoms

    return NativeData(atoms=atoms, relax=relax, spins=spins)





def read_lat_array(lines):
    lat_array = []
    for lat_vec in lines[:3]:  # lattice vectors in next 3 lines
        vec = np.array([float(a) for a in lat_vec.split()])
        lat_array.append(vec / np.linalg.norm(vec))  # normalize
    return np.array(lat_array)


def read_inpt_cell(path):
    lat_array = []
    cell = []
    with open(path, "r") as inptfile:
        contents = inptfile.read()

    lines = contents.splitlines()
    for i, line in enumerate(lines):
        if "CELL" in line:
            cell = np.array([float(a) for a in line.split()[1:]])
        elif "LATVEC" in line:
            lat_array = read_lat_array(lines[i + 1 :])

    if len(cell) == 0:
        raise Exception(f"no cell found in {path}")
    if len(lat_array) == 0:
        lat_array = np.eye(3)

    return (lat_array.T * cell).T * Bohr


def validate_atom_block(block):
    assert "ATOM_TYPE" in block and block["ATOM_TYPE"]
    assert ("COORD" in block) ^ ("COORD_FRAC" in block)  # XOR
    assert "N_TYPE_ATOM" in block
    natoms = int(block["N_TYPE_ATOM"])
    for key in ["COORD", "COORD_FRAC", "RELAX", "SPIN"]:
        assert key not in block or natoms == len(block[key])


def reorder(original, order):
    res = original.copy()
    for oldi, newi in enumerate(order):
        res[newi] = original[oldi]
    return res


def constraints_from_relax(constraints):
    """
    This is just a helper function to translate from a list of
    constraints from SPARC to a set of ASE constraints

    Parameters:
        constraints (list):
            a list for lists, either empty or in the form [x,y,z]
            i.e. [[],[0,0,1],[1,1,1]]

    returns:
        cons_list (list)
            a list of ase contraint classes
    """
    cons_list = []
    fix_atoms = []
    for i, constraint in enumerate(constraints):
        if constraint.count(1) == 3 or constraint == []:
            continue
        elif constraint.count(1) == 0:
            fix_atoms.append(i)
        elif constraint.count(1) == 1:
            cons_list.append(FixedLine(i, constraint))
        elif constraint.count(1) == 2:
            tmp = [int(not a) for a in constraint]  # flip ones and zeros
            cons_list.append(FixedPlane(i, tmp))
    if len(fix_atoms) != 0:
        cons_list.append(FixAtoms(fix_atoms))
    return cons_list


def format_float(num):
    # the space in front of the decimal place replaces the negative with a space
    # in positive numbers
    return f"{num: .15f}"


def format_cell(cell):
    key = "#CELL:"
    comps = [format_float(np.linalg.norm(comp) / Bohr) for comp in cell]
    return "  ".join([key] + comps)


def format_latvec(cell):
    def format_row(comp):
        comp_n = [format_float(c) for c in comp / np.linalg.norm(comp)]
        return " ".join(["#"] + comp_n)

    key = "#LATVEC"
    rows = [format_row(r) for r in cell]
    return "\n".join([key] + rows)


def format_pbc(pbc):
    key = "#PBC:"
    return " ".join([key] + [str(c) for c in pbc])


def relax_from_constraint(constraint):
    """returns dict of {atom_index: relax_dimensions} for the given constraint"""
    if isinstance(constraint, FixAtoms):
        dimensions = [0] * 3
        expected_free = 0
    elif isinstance(constraint, (FixedLine, FixedPlane)):
        dimensions = [int(bool(a)) for a in constraint.dir]
        expected_free = 1
    elif isinstance(constraint, FixedPlane):
        dimensions = [int(not a) for a in constraint.dir]
        expected_free = 2
    else:
        warnings.warn(
            "The constraint type {} is not supported by"
            " the .ion format. This constraint will be"
            " ignored"
        )
        return {}
    if dimensions.count(1) != expected_free:
        warnings.warn(
            "The .ion filetype can only support freezing entire "
            f"dimensions (x,y,z). The {type(constraint).__name__} "
            "constraint on this atoms Object is being converted to "
            "allow movement in any directions in which the atom had "
            "some freedom to move"
        )
    return {i: dimensions for i in constraint.get_indices()}  # atom indices


def relax_from_all_constraints(constraints, natoms):
    """converts ASE atom constraints to SPARC relaxed dimensions for the atoms"""
    relax = [[1] * 3] * natoms  # assume relaxed in all dimensions for all atoms
    for c in constraints:
        for atom_index, rdims in relax_from_constraint(c).items():
            # 'and' current relaxed dimensions with relaxed dimensions from
            # constraint c (constrain if either is constrained)
            relax[atom_index] = [int(o and n) for o, n in zip(relax[i], rdims)]
    return relax


def format_atom_block(atoms: List[Atom], pseudo_path: str, relax=[], scaled=True):
    # All atoms here should be same type
    assert len({atom.symbol for atom in atoms}) == 1
    atom = atoms[0]
    res = ""
    res += f"ATOM_TYPE: {atom.symbol}\n"
    res += f"N_TYPE_ATOM: {len(atoms)}\n"
    res += f"PSEUDO_POT: {pseudo_path}\n"
    res += f"ATOMIC_MASS: {atom.mass}\n"

    def coord_row(coord, index) -> str:
        nums = [format_float(x) for x in coord]
        coord_string = "    ".join(nums)
        return f"    {coord_string}   # index {index}"

    def relax_row(i):
        return " ".join(str(a) for a in relax[i])

    if scaled:
        # for some reason atom.scaled_position are rounded differently from
        # the scaled positions below. Using this version for total backwards
        # compatability
        sp = atoms[0].atoms.get_scaled_positions(wrap=False)
        res += "COORD_FRAC:\n"
        res += "\n".join(coord_row(sp[a.index], a.index) for a in atoms)
    else:
        res += "COORD:\n"
        res += "\n".join(coord_row(a.position / Bohr, a.index) for a in atoms)

    if len(relax) != 0:
        res += "\nRELAX:\n"
        res += "\n".join(relax_row(a.index) for a in atoms)

    magmoms = [a.magmom for a in atoms]
    if any(magmoms):
        res += "\nSPIN:\n"
        res += "\n".join(format_float(m) for m in magmoms)
    return res


def find_pseudo_path(element: str, pseudo_dir):
    suffix = element + ".pot"
    filenames = [a for a in os.listdir(pseudo_dir) if a.endswith(suffix)]
    if len(filenames) == 0:
        raise Exception("File not found")

    filename = filenames[0]
    if len(filenames) > 1:
        warnings.warn(
            "Multiple psudopotentials detected for "
            f"{element} ({filenames}). Using {filename}."
        )
    return os.path.join(pseudo_dir, filename)


def copy_psp(srcpath: str, dstdir: str):
    filename = os.path.basename(srcpath)
    dstpath = os.path.join(dstdir, filename)
    if srcpath != dstpath:
        shutil.copyfile(srcpath, dstpath)
