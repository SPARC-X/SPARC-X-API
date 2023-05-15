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

from .utils import read_block_input, bisect_and_strip

from ..api import SparcAPI
import re
from datetime import datetime

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
    print(sparc_version)
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
    date_version = datetime.strptime(date_str, "%b %d %Y").strftime("%Y.%m.%d")
    return date_version


def _read_input_params(contents, validator=defaultAPI):
    """Parse the Input parameters and Paral"""
    lines = "\n".join(
        _get_block_text(contents, "Input parameters")
        + _get_block_text(contents, "Parallelization")
    ).split("\n")
    print(lines)
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
    convergence_info = _get_block_text(
        contents, r"Self Consistent Field \(SCF.*?\)"
    )
    results_info = _get_block_text(contents, "Energy and force calculation")

    if len(convergence_info) != len(results_info):
        # TODO: change to another exception name
        raise ValueError(
            "Error, length of convergence information and energy calculation are different!"
        )
    n_steps = len(convergence_info)
    steps = []
    for i, step in enumerate(zip(convergence_info, results_info)):
        current_step = {"scf_step": i}
        conv, res = step
        # TODO: add support for convergence fields
        conv_lines = conv.splitlines()
        conv_header = re.split(r"\s{3,}", conv_lines[0])
        # omit the last line which is just a checker
        conv_array = np.genfromtxt(conv_lines[1:-1], dtype=float)
        # TODO: the meaning of the header should me split to the width

        conv_dict = {}
        for i, field in enumerate(conv_header):
            field = field.split("(")[0].strip().lower()
            value = conv_array[:, i]
            if "free energy" in field:
                value *= Hartree
            conv_dict[field] = value

        current_step["convergence"] = conv_dict

        {"header": conv_header, "values": conv_array}

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
    raise NotImplementedError(
        "Writing output file from python-api not supported!"
    )
