"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)

This file has been heavily modified since SPARC 0.1

TODO: more descriptions about this file io parser
"""
import re
from datetime import datetime
from warnings import warn

import numpy as np
from ase.units import Bohr, GPa, Hartree

# Safe wrappers for both string and fd
from ase.utils import reader, writer

from ..api import SparcAPI
from .utils import bisect_and_strip, read_block_input

# TODO: should allow user to select the api
defaultAPI = SparcAPI()


@reader
def _read_out(fileobj):
    """
    Read the .out file content

    The output file is just stdout. The blocks are read using re-patterns rather than the way .static / .geopt or .aimd are parsed
    """
    contents = fileobj.read()
    sparc_version = _read_sparc_version(contents[:4096])
    # print(sparc_version)
    # TODO: use the sparc version to construct the API
    output_dict = {"sparc_version": sparc_version}
    # Combine the input parameters and parallelization sections
    output_dict["parameters"] = _read_input_params(contents)

    # Parse the Initialization and timing info, and if calculation
    # successfully finished
    # Note: not all information are converted!
    output_dict["run_info"] = _read_run_info(contents)
    # List of scf information,
    # including scf convergence, energy etc
    output_dict["ionic_steps"] = _read_scfs(contents)
    return {"out": output_dict}


def _read_sparc_version(header):
    """Read the sparc version from the output file header.

    This function should live outside the _read_output since some other functions may use it

    TODO: combine it with the version from initialization.c
    """
    pattern_version = r"SPARC\s+\(\s*?version(.*?)\)"
    match = re.findall(pattern_version, header)
    if len(match) != 1:
        warn("Header does not contain SPARC version information!")
        return None
    date_str = match[0].strip().replace(",", " ")
    # Accept both abbreviate and full month name
    try:
        date_version = datetime.strptime(date_str, "%B %d %Y").strftime("%Y.%m.%d")
    except ValueError:
        try:
            date_version = datetime.strptime(date_str, "%b %d %Y").strftime("%Y.%m.%d")
        except ValueError:
            warn("Cannot fetch SPARC version information!")
            date_version = None
    return date_version


def _read_input_params(contents, validator=defaultAPI):
    """Parse the Input parameters and Paral"""
    lines = "\n".join(
        _get_block_text(contents, "Input parameters")
        + _get_block_text(contents, "Parallelization")
    ).split("\n")
    # print(lines)
    params = read_block_input(lines, validator=validator)
    return params


def _read_run_info(contents):
    """Parse the run info sections
    Note due to the complexity of the run info,
    the types are not directly converted
    """
    lines = "\n".join(
        _get_block_text(contents, "Timing info")
        + _get_block_text(contents, "Initialization")
    ).split("\n")
    block_dict = {"raw_info": lines}
    # Select key fields to store
    for line in lines:
        if ":" not in line:
            continue
        key, value = bisect_and_strip(line, ":")
        key = key.lower()
        if key in block_dict:
            if key not in ("pseudopotential",):
                warn(
                    f"Key {key} from run information appears multiple times in your outputfile!"
                )
            # For keys like pseudopotential, we make it a list
            else:
                origin_value = list(block_dict[key])
                value = origin_value + [value]

        block_dict[key] = value
    return block_dict


def _read_scfs(contents):
    """Parse the ionic steps

    Return:
    List of ionic steps information


    """
    convergence_info = _get_block_text(contents, r"Self Consistent Field \(SCF.*?\)")
    results_info = _get_block_text(contents, "Energy and force calculation")

    # Should not happen
    if len(convergence_info) > len(results_info) + 1:
        raise ValueError(
            "Error, length of convergence information and energy calculation mismatch!"
        )
    elif len(convergence_info) == len(results_info) + 1:
        warn("Last ionic SCF has not finished! The results may be incomplete")
    else:
        pass

    # Stick to the convergence information as the main section
    n_steps = len(convergence_info)
    steps = []
    # for i, step in enumerate(zip(convergence_info, results_info)):
    for i in range(n_steps):
        conv = convergence_info[i]
        # Solution for incomplete calculations
        if i >= len(results_info):
            res = ""  # Empty lines
        else:
            res = results_info[i]
        current_step = {"scf_step": i}
        # TODO: add support for convergence fields
        conv_lines = conv.splitlines()
        # conv_header is normally 4-column table
        conv_header = re.split(r"\s{3,}", conv_lines[0])

        scf_sub_steps = []
        # For ground-state calculations, the output will be only 1 block
        # For hybrid (HSE/PBE0) calculations the EXX loops will also be included
        # General rule: we search for the line "Total number of SCF: N", read back N(+1) lines
        for lino, line in enumerate(conv_lines):
            if "Total number of SCF:" not in line:
                continue
            scf_num = int(line.split(":")[-1])
            conv_array = np.genfromtxt(
                [
                    l
                    for l in conv_lines[lino - scf_num : lino]
                    if l.split()[0].isdigit()
                ],
                dtype=float,
                ndmin=2,
            )
            conv_dict = {}
            for i, field in enumerate(conv_header):
                field = field.strip()
                value = conv_array[:, i]
                # TODO: re-use the value conversion function in res part
                if "Ha/atom" in field:
                    value *= Hartree
                    field.replace("Ha/atom", "eV/atom")
                if "Iteration" in field:
                    value = value.astype(int)
                conv_dict[field] = value
            # Determine if the current block is a ground-state or EXX
            name_line = conv_lines[lino - scf_num - 1]
            if "Iteration" in name_line:
                name = "ground state"
            else:
                name = name_line

            conv_dict["name"] = name
            scf_sub_steps.append(conv_dict)

        current_step["convergence"] = scf_sub_steps

        res = res.splitlines()
        for line in res:
            if ":" not in line:
                continue
            key, value = bisect_and_strip(line, ":")
            key = key.lower()
            if key in current_step:
                warn(
                    f"Key {key} appears multiples in one energy / force calculation, your output file may be incorrect."
                )
            # Conversion of values are relatively easy
            pattern_value = r"([+\-\d.Ee]+)\s+\((.*?)\)"
            match = re.findall(pattern_value, value)
            raw_value, unit = float(match[0][0]), match[0][1]
            if unit == "Ha":
                converted_value = raw_value * Hartree
                converted_unit = "eV"
            elif unit == "Ha/atom":
                converted_value = raw_value * Hartree
                converted_unit = "eV/atom"
            elif unit == "Ha/Bohr":
                converted_value = raw_value * Hartree / Bohr
                converted_unit = "eV/Angstrom"
            elif unit == "GPa":
                converted_value = raw_value * GPa
                converted_unit = "eV/Angstrom^3"
            elif unit == "sec":
                converted_value = raw_value * 1
                converted_unit = "sec"
            elif unit == "Bohr magneton":
                converted_value = raw_value
                converted_unit = "Bohr magneton"
            else:
                warn(f"Conversion for unit {unit} unknown! Treat as unit")
                converted_value = raw_value
                converted_unit = unit
            current_step[key] = {
                "value": converted_value,
                "unit": converted_unit,
            }
        steps.append(current_step)
    return steps


def _get_block_text(text, block_name):
    """Get an output 'block' with a specific block name

    the outputs are not line-split
    """
    # Add the ending separator so matching is possible for partial-complete
    # .out file from socket calculations
    text = text + ("=" * 68) + "\n"
    pattern_block = (
        r"[\*=]{50,}\s*?\n\s*?BLOCK_NAME\s*?\n[\*=]{50,}\s*\n(.*?)[\*=]{50,}"
    )
    pattern = pattern_block.replace("BLOCK_NAME", block_name)
    match = re.findall(pattern, text, re.DOTALL | re.MULTILINE)
    if len(match) == 0:
        warn(f"Block {block_name} cannot be parsed from current text!")
    return match


@writer
def _write_out(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing output file from SPARC-X-API not supported!")
