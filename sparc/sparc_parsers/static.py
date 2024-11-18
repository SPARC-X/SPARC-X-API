"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)

This file has been heavily modified since SPARC 0.1

TODO: more descriptions about this file io parser
"""
from warnings import warn

import numpy as np
from ase.units import Bohr, GPa, Hartree

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from ..api import SparcAPI
from .utils import strip_comments

# TODO: should allow user to select the api
# defaultAPI = SparcAPI()


@reader
def _read_static(fileobj):
    """
    Read the .static file content

    Each .static file should only host 1 or more images
    (if socket mode is enabled), but the output may vary
    a lot depending on the flags (e.g. PRINT_ATOMS, PRINT_FORCES etc)
    """
    contents = fileobj.read()
    data, comments = strip_comments(contents)
    # Most static files should containe the Atom positions lines
    # this can happen for both a single- or multi-image file
    step_bounds = [i for i, x in enumerate(data) if "Atom positions" in x] + [len(data)]
    raw_static_steps = [
        data[start:end] for start, end in zip(step_bounds[:-1], step_bounds[1:])
    ]
    # In some cases (e.g. PRINT_ATOMS=0), the static file may not contain Atom positions
    # All existing lines will be regarded as one step
    if len(raw_static_steps) == 0:
        raw_static_steps = [data]
    static_steps = [_read_static_step(step) for step in raw_static_steps]
    return {"static": static_steps}


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
    elif "Net magnetization" in header_name:
        name = "net_magnetization"
    elif "Atomic magnetization" in header_name:
        name = "atomic_magnetization"
    elif "Stress (GPa)" in header_name:
        name = "stress"
    elif "Stress equiv." in header_name:
        name = "stress_equiv"
    elif "Stress (Ha/Bohr)" in header_name:
        name = "stress_1d"
    elif "Stress (Ha/Bohr**2)" in header_name:
        name = "stress_2d"
    elif "Fractional coordinates" in header_name:
        # Fractional coordinates of Si -- > name=coord_frac symbol="Si"
        name = "coord_frac"
        symbol = header_name.split("of")[1].strip()
        clean_array = value.reshape((-1, 3))
        value = {"value": clean_array, "symbol": symbol}
    # Exclusive to the socket mode
    elif "Lattice (Bohr)" in header_name:
        name = "lattice"
        value = value.reshape((3, 3))
    else:
        name = header_name.strip()

    return {"name": name, "value": value}


def _read_static_step(step):
    """Parse all the lines in one step and compose a dict containing sanitized blocks
    Args:
        step (list): Lines of raw lines in one step
    """
    separator = "*" * 60  # Make the separator long enough
    # Clean up boundary lines
    data = [
        line
        for line in step
        if ("Atom positions" not in line) and (separator not in line)
    ]
    block_bounds = [i for i, x in enumerate(data) if ":" in x] + [len(data)]
    raw_blocks = [
        data[start:end] for start, end in zip(block_bounds[:-1], block_bounds[1:])
    ]
    # import pdb; pdb.set_trace()
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
                # import pdb; pdb.set_trace()
                coord_dict["coord_frac"] = np.vstack(
                    [coord_dict["coord_frac"], coord_f]
                )

        elif name == "free energy":
            value = raw_value * Hartree
        elif name == "forces":
            value = raw_value * Hartree / Bohr
        elif name == "atomic_magnetization":
            value = raw_value
        elif name == "net_magnetization":
            value = raw_value
        elif name == "stress":
            # Stress is in eV/Ang^3, may need to convert to Virial later when cell is known
            # For low-dimension stress info, use stress_equiv
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
        elif name == "stress_equiv":
            # Only store the size up to the max. periodic directions,
            # let the atom parser decide how to transform the matrix
            value = raw_value * GPa
        elif name == "stress_1d":
            value = raw_value * Hartree / Bohr
        elif name == "stress_2d":
            value = raw_value * Hartree / (Bohr**2)
        elif name == "lattice":
            value = raw_value * Bohr

        # Non-frac coord
        if value is not None:
            block_dict[name] = value
    # Finally, update the atomic positions
    # TODO: should we keep a default?
    if coord_dict != {}:
        block_dict["atoms"] = coord_dict

    return block_dict


def _add_cell_info(static_steps, cell=None):
    """Use the cell information to convert positions
    if lattice exists in each step, use it to convert coord_frac
    else use the external cell (for example from inpt file)
    Args:
        static_steps: raw list of steps
        cell: external lattice information
    """
    new_steps = []
    for step in static_steps:
        new_step = step.copy()
        if "lattice" in step:
            lat = step["lattice"]
        elif cell is not None:
            lat = cell
        else:
            lat = None

        if (lat is not None) and (step.get("atoms", None) is not None):
            coord_frac = new_step["atoms"]["coord_frac"]
            coord_cart = np.dot(coord_frac, lat)
            new_step["atoms"]["coord"] = coord_cart
        new_steps.append(new_step)
    return new_steps


@writer
def _write_static(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing static file from SPARC-X-API not supported!")
