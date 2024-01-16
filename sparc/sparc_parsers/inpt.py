import numpy as np
from ase.units import Bohr

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from ..api import SparcAPI
from .utils import read_block_input, strip_comments

defaultAPI = SparcAPI()


@reader
def _read_inpt(fileobj, validator=defaultAPI):
    contents = fileobj.read()
    # label = get_label(fileobj, ".ion")
    data, comments = strip_comments(contents)
    # We do not read the cell at this time!

    # find the index for all atom type lines. They should be at the
    # top of their block
    inpt_blocks = read_block_input(data, validator=validator)
    return {"inpt": {"params": inpt_blocks, "comments": comments}}


@writer
def _write_inpt(fileobj, data_dict, validator=defaultAPI):
    if "inpt" not in data_dict:
        raise ValueError("Your dict does not contain inpt section!")

    inpt_dict = data_dict["inpt"]

    if "params" not in inpt_dict:
        raise ValueError("Input dict for inpt file does not have `params` field!")

    comments = inpt_dict.get("comments", [])
    banner = "Input File Generated By SPARC ASE Calculator"
    if len(comments) == 0:
        comments = [banner]
    elif "ASE" not in comments[0]:
        comments = [banner] + comments
    for line in comments:
        fileobj.write(f"# {line}\n")
    fileobj.write("\n")
    params = inpt_dict["params"]
    for key, val in params.items():
        # TODO: can we add a multiline argument?
        val_string = validator.convert_value_to_string(key, val)
        if (val_string.count("\n") > 0) or (
            key
            in [
                "LATVEC",
            ]
        ):
            output = f"{key}:\n{val_string}\n"
        else:
            output = f"{key}: {val_string}\n"
        fileobj.write(output)
    return


def _inpt_cell_to_ase_cell(data_dict):
    """Convert the inpt cell convention to a real cell (in ASE Angstrom unit)

    Arguments:
    inpt_blocks: an already validated inpt file blocks dict
                 (i.e. parsed by _read_inpt)

    Returns:
    cell in ASE-unit
    """
    inpt_blocks = data_dict["inpt"]["params"]
    if ("CELL" in inpt_blocks) and ("LATVEC_SCALE" in inpt_blocks):
        # TODO: customize the exception class
        # TODO: how do we convert the rule from doc?
        raise ValueError("LATVEC_SCALE and CELL cannot be specified simultaneously!")

    # if "CELL" in inpt_blocks:
    #     cell = np.eye(inpt_blocks["CELL"]) * Bohr
    if "LATVEC" not in inpt_blocks:
        if ("CELL" in inpt_blocks) or ("LATVEC_SCALE" in inpt_blocks):
            lat_array = np.eye(3) * Bohr
        else:
            raise KeyError(
                "LATVEC is missing in inpt file and no CELL / LATVEC_SCALE provided!"
            )
    else:
        lat_array = np.array(inpt_blocks["LATVEC"]) * Bohr

    # LATVEC_SCALE: just multiplies
    if "LATVEC_SCALE" in inpt_blocks:
        scale = np.array(inpt_blocks["LATVEC_SCALE"])
        cell = (lat_array.T * scale).T

    # CELL: the lengths are in the LATVEC directions
    # TODO: the documentation about CELL is a bit messy. Is CELL always orthogonal?
    # Anyway the lat_array when CELL is none should be ok
    elif "CELL" in inpt_blocks:
        scale = np.array(inpt_blocks["CELL"])
        unit_lat_array = (
            lat_array / np.linalg.norm(lat_array, axis=1, keepdims=True) * Bohr
        )
        cell = (unit_lat_array.T * scale).T
    else:
        cell = lat_array
    return cell
