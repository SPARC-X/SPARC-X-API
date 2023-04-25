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
from ase.units import Bohr, Hartree, eV, GPa
from ase.constraints import FixAtoms, FixedLine, FixedPlane


# Safe wrappers for both string and fd
from ase.utils import reader, writer

from .utils import (
    get_label,
    strip_comments,
    bisect_and_strip,
    make_reverse_mapping,
)

from ..inputs import SparcInputs
import textwrap
import re
from datetime import datetime

# TODO: should allow user to select the api
defaultAPI = SparcInputs()


@reader
def _read_out(fileobj):
    """
    Read the .out file content

    The output file is just a formatted stdout, so there are many chances that 

    Each .static file should only host 1 image (as least per now), but the output may vary
    a lot depending on the flags (e.g. PRINT_ATOMS, PRINT_FORCES etc)
    """
    contents = fileobj.read()
    sparc_version = _read_sparc_version(contents[:4096])
    print(sparc_version)


def _read_sparc_version(header):
    """Read the sparc version from the output file header.

    This function should live outside the _read_output since some other functions may use it

    TODO: combine it with the version from initialization.c
    """
    pattern_version = r"SPARC\s+\(\s*?version(.*?)\)"
    match = re.findall(pattern_version, header)
    if len(match) != 1:
        warn(
            "Header does not contain SPARC version information!"
        )
        return None
    date_str = match[0].strip().replace(",", " ")
    date_version = datetime.strptime(date_str, "%b %d %Y").strftime("%Y.%m.%d")
    return date_version
    
    

@writer
def _write_out(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing output file from python-api not supported!")
