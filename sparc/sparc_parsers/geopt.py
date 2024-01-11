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
defaultAPI = SparcAPI()


@reader
def _read_geopt(fileobj):
    """
    Parse the geopt information
    Each geopt is similar to the static block, except that the field name is started by
    ':'
    The relaxations are separated by ':RELAXSTEP:' seperators
    """
    contents = fileobj.read()
    # label = get_label(fileobj, ".ion")
    # The geopt comments are simply discarded
    data, comments = strip_comments(contents)

    # find the index for all atom type lines. They should be at the
    # top of their block
    step_bounds = [i for i, x in enumerate(data) if ":RELAXSTEP:" in x] + [len(data)]
    raw_geopt_blocks = [
        data[start:end] for start, end in zip(step_bounds[:-1], step_bounds[1:])
    ]
    geopt_steps = [_read_geopt_step(step) for step in raw_geopt_blocks]

    return {"geopt": geopt_steps}


def _read_geopt_step(raw_step_text):
    """Parse a geopt step and compose the data dict

    Arguments
    raw_step_text: list of lines within the step
    """
    header, body = raw_step_text[0], raw_step_text[1:]
    if ":RELAXSTEP:" not in header:
        raise ValueError("Wrong geopt format! The :RELAXSTEP: label is missing.")
    # Geopt file uses 1-indexed step names, convert to 0-indexed
    step = int(header.split(":RELAXSTEP:")[-1]) - 1
    bounds = [i for i, x in enumerate(body) if ":" in x] + [len(body)]
    blocks = [body[start:end] for start, end in zip(bounds[:-1], bounds[1:])]
    data = {}
    for block in blocks:
        header_block, body_block = block[0], block[1:]
        header_name = header_block.split(":")[1]
        header_data = header_block.split(":")[-1].strip()
        if len(header_data) > 0:
            block_raw_data = [header_data] + body_block
        else:
            block_raw_data = body_block
        # import pdb; pdb.set_trace()
        raw_value = np.genfromtxt(block_raw_data, dtype=float)
        if "R(Bohr)" in header_name:
            name = "positions"
            value = raw_value.reshape((-1, 3)) * Bohr
        elif "E(Ha)" in header_name:
            name = "energy"
            value = float(raw_value) * Hartree
        elif "F(Ha/Bohr)" in header_name:
            name = "forces"
            value = raw_value.reshape((-1, 3)) * Hartree / Bohr
        elif "CELL" in header_name:
            name = "cell"
            value = raw_value * Bohr
        elif "VOLUME" in header_name:
            name = "volume"
            value = raw_value * Bohr**3
        elif "LATVEC" in header_name:
            # TODO: the LATVEC is ambiguous. Are the results a unit cell, or full cell?
            name = "latvec"
            value = raw_value * Bohr
        elif "STRESS" in header_name:
            # Stress handling in geopt output can be different
            # on low-dimensional systems. If the stress matrix is 3x3,
            # the unit is GPa, while lower dimensional stress matrices
            # are using Hartree / Bohr**2 or Hartree / Bohr
            dim = raw_value.shape[0]
            if dim == 3:
                name = "stress"
                stress_ev_a3 = raw_value * GPa
                # Standard stress value, use Voigt representation
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
            elif dim == 2:
                name = "stress_2d"
                value = raw_value * Hartree / Bohr**2
            elif dim == 1:
                name = "stress_1d"
                value = raw_value * Hartree / Bohr
            else:
                raise ValueError("Incorrect stress matrix dimension!")
        else:
            warn(
                f"Field {header_name} is not known to geopt! I'll use the results as is."
            )
            name = header_name
            value = raw_value
        data[name] = value
    # Special treatment for latvec & cell
    if ("cell" in data) and ("latvec" in data):
        # TODO: check
        cell_, lat_ = data["cell"], data["latvec"]
        unit_lat = lat_ / np.linalg.norm(lat_, axis=1, keepdims=True)
        cell = (unit_lat.T * cell_).T
        data["ase_cell"] = cell
    data["step"] = step
    return data


@writer
def _write_geopt(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing geopt file from SPARC-X-API not supported!")
