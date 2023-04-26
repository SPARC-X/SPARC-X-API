"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)

This file has been heavily modified since SPARC 0.1

TODO: more descriptions about this file io parser
"""
from warnings import warn

import numpy as np
from ase.units import Bohr, Hartree, GPa

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from .utils import (
    strip_comments,
)

from ..inputs import SparcInputs

# TODO: should allow user to select the api
defaultAPI = SparcInputs()


@reader
def _read_static(fileobj):
    """
    Read the .static file content

    Each .static file should only host 1 image (as least per now), but the output may vary
    a lot depending on the flags (e.g. PRINT_ATOMS, PRINT_FORCES etc)
    """
    contents = fileobj.read()
    data, comments = strip_comments(contents)
    # Like .ion file, but split by the separator ":"
    block_bounds = [i for i, x in enumerate(data) if ":" in x] + [len(data)]
    # blocks = [read_static_block(data[start:end]) for start, end in zip(block_bounds[:-1], block_bounds[1:])]
    raw_blocks = [
        data[start:end] for start, end in zip(block_bounds[:-1], block_bounds[1:])
    ]
    static_dict = read_static_blocks(raw_blocks)
    return {"static": static_dict}


def _read_static_block(raw_block):
    """Parse ONE static data block, this will return dict with keys:
    {"name": PARAM, "value": value}

    Arguments:
    raw_block: un-parsed block as list of strings
    """
    header, body = raw_block[0], raw_block[1:]
    header_name, header_rest = header.split(":")
    if len(header_rest.strip()) > 0:
        body = [header_rest.strip()] + body

    try:
        value = np.genfromtxt(body, dtype=float)
        if np.isnan(value).any():
            warn(
                (
                    f"Field contains data that are not parsable by numpy!\n"
                    f"Contents are: {body}"
                )
            )
    except Exception:
        value = body

    # name = None
    if "Total free energy" in header_name:
        name = "free energy"
    elif "Atomic forces" in header_name:
        name = "forces"
    elif "Stress" in header_name:
        name = "stress"
    elif "Fractional coordinates" in header_name:
        # Fractional coordinates of Si -- > name=coord_frac symbol="Si"
        name = "coord_frac"
        symbol = header_name.split("of")[1].strip()
        clean_array = value.reshape((-1, 3))
        value = {"value": clean_array, "symbol": symbol}
    else:
        name = header_name.strip()

    return {"name": name, "value": value}


def read_static_blocks(raw_blocks):
    """Read all blocks from the static file and compose a dict"""
    block_dict = {}
    coord_dict = {}
    block_contents = [_read_static_block(block) for block in raw_blocks]
    for bc in block_contents:
        name, raw_value = bc["name"], bc["value"]
        # Coord frac needs to be collected in all positions
        if name == "coord_frac":
            value = None
            symbol, coord_f = raw_value["symbol"], raw_value["value"]
            pos_count = len(coord_f)
            if coord_dict == {}:
                coord_dict.update(
                    {
                        "symbols": [
                            symbol,
                        ]
                        * pos_count,
                        "coord_frac": coord_f,
                    }
                )
            else:
                coord_dict["symbols"] += [
                    symbol,
                ] * pos_count
                coord_dict["coord_frac"] += np.vstack(
                    [coord_dict["coord_frac"], coord_f]
                )

        elif name == "free energy":
            value = raw_value * Hartree
        elif name == "forces":
            value = raw_value * Hartree / Bohr
        elif name == "stress":
            # Stress is in eV/Ang^3, may need to convert to Virial later when cell is known
            stress_ev_a3 = raw_value * GPa
            if stress_ev_a3.shape != (3, 3):
                raise ValueError("Stress from static file is not a 3x3 matrix!")
            # make the stress in voigt notation
            # TODO: check the order!
            value = np.array(
                [
                    stress_ev_a3[0, 0],
                    stress_ev_a3[1, 1],
                    stress_ev_a3[2, 2],
                    stress_ev_a3[1, 2],
                    stress_ev_a3[0, 2],
                    stress_ev_a3[0, 1],
                ]
            )

        # Non-frac coord
        if value is not None:
            block_dict[name] = value
    # Finally, update the atomic positions
    # TODO: should we keep a default?
    if coord_dict != {}:
        block_dict["atoms"] = coord_dict

    return block_dict


def _add_cell_info(static_dict, cell=None):
    """When cell information is available, convert

    The cell should already be in angstrom

    1) the cell position from fractional to cartesian
    2) Anything else?
    """
    new_dict = static_dict.copy()
    block_dict = new_dict["static"]
    if (cell is None) or (block_dict.get("atoms", None) is None):
        return new_dict
    coord_frac = block_dict["atoms"]["coord_frac"]
    coord_cart = np.dot(coord_frac, cell)
    new_dict["static"]["atoms"]["coord"] = coord_cart
    return new_dict


@writer
def _write_static(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing static file from python-api not supported!")
