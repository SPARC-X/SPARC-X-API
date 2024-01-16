"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)

This file has been heavily modified since SPARC 0.1

TODO: more descriptions about this file io parser
"""
from warnings import warn

import numpy as np
from ase.units import AUT, Angstrom, Bohr, GPa, Hartree, fs

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from ..api import SparcAPI
from .utils import strip_comments


@reader
def _read_aimd(fileobj):
    """Parse the aimd information Each geopt is similar to the static
    block, except that the field name is started by ':' The
    relaxations are separated by ':MDSTEP:' seperators

    """
    contents = fileobj.read()
    # label = get_label(fileobj, ".ion")
    # The geopt comments are simply discarded
    stripped, comments = strip_comments(contents)
    # Do not include the description lines
    data = [line for line in stripped if ":Desc" not in line]

    # find the index for all atom type lines. They should be at the
    # top of their block
    step_bounds = [i for i, x in enumerate(data) if ":MDSTEP:" in x] + [len(data)]
    raw_aimd_blocks = [
        data[start:end] for start, end in zip(step_bounds[:-1], step_bounds[1:])
    ]
    aimd_steps = [_read_aimd_step(step) for step in raw_aimd_blocks]

    return {"aimd": aimd_steps}


def _read_aimd_step(raw_aimd_text):
    """Parse a geopt step and compose the data dict

    Arguments
    raw_aimd_text: list of lines within the step

    Most values are just presented in their output format,
    higher level function calling _read_aimd_step and _read_aimd
    should implement how to use the values,
    e.g. E_tot = E_tot_per_atom * N_atoms

    """
    header, body = raw_aimd_text[0], raw_aimd_text[1:]
    if ":MDSTEP:" not in header:
        raise ValueError("Wrong aimd format! The :MDSTEP: label is missing.")
    # Geopt file uses 1-indexed step names, convert to 0-indexed
    step = int(header.split(":MDSTEP:")[-1]) - 1
    print("Step ", step)
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
        # The type definitions from MD may be treated from API again?
        if header_name == "R":
            name = "positions"
            value = raw_value.reshape((-1, 3)) * Bohr
        elif header_name == "V":
            name = "velocities"
            value = raw_value.reshape((-1, 3)) * Bohr / AUT / (Angstrom / fs)
        elif header_name == "F":
            name = "forces"
            value = raw_value.reshape((-1, 3)) * Hartree / Bohr
        elif header_name == "MDTM":
            # This is not the md integration time!
            name = "md_walltime"
            value = float(raw_value)
        elif header_name == "TEL":
            name = "electron temp"
            value = float(raw_value)
        elif header_name == "TIO":
            name = "ion temp"
            value = float(raw_value)
        elif header_name == "TEN":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "total energy per atom"
            value = float(raw_value) * Hartree
        elif header_name == "KEN":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "kinetic energy per atom"
            value = float(raw_value) * Hartree
        elif header_name == "KENIG":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "kinetic energy (ideal gas) per atom"
            value = float(raw_value) * Hartree
        elif header_name == "FEN":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "free energy per atom"
            value = float(raw_value) * Hartree
        elif header_name == "UEN":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "internal energy per atom"
            value = float(raw_value) * Hartree
        elif header_name == "TSEN":
            # Note it's the total energy per atom!
            # TODO: shall we convert to ase fashion?
            name = "entropy*T per atom"
            value = float(raw_value) * Hartree
        elif header_name == "STRESS":
            # Same rule as STRESS in geopt
            # no conversion to Voigt form yet
            dim = raw_value.shape[0]
            if dim == 3:
                name = "stress"
                value = raw_value * GPa
            elif dim == 2:
                name = "stress_2d"
                value = raw_value * Hartree / Bohr**2
            elif dim == 1:
                name = "stress_1d"
                value = raw_value * Hartree / Bohr
            else:
                raise ValueError("Incorrect stress matrix dimension!")
        elif header_name == "STRIO":
            # Don't do the volume conversion now
            name = "stress (ion-kinetic)"
            value = raw_value * GPa
        elif header_name == "PRES":
            # Don't do the volume conversion now
            name = "pressure"
            value = raw_value * GPa
        elif header_name == "PRESIO":
            # Don't do the volume conversion now
            name = "pressure (ion-kinetic)"
            value = raw_value * GPa
        elif header_name == "PRESIG":
            # Don't do the volume conversion now
            name = "pressure (ideal gas)"
            value = raw_value * GPa
        elif header_name in ("AVGV", "MAXV", "MIND"):
            warn(f"MD output keyword {header_name} will not be parsed.")
            value = None
        else:
            warn(f"MD output keyword {header_name} not known to SPARC. " "Ignore.")
            value = None
        if value is not None:
            data[name] = value
    data["step"] = step
    return data


@writer
def _write_aimd(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing aimd file from SPARC-X-API " "not supported!")
