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

from .utils import (
    get_label,
    strip_comments,
    bisect_and_strip,
    read_block_input,
    make_reverse_mapping,
)

from ..inputs import SparcInputs
import textwrap

# TODO: should allow user to select the api
defaultAPI = SparcInputs()


@reader
def _read_aimd(fileobj):
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
    sort, resort, new_comments = _read_sort_comment(comments)

    # find the index for all atom type lines. They should be at the top of their block
    atom_type_bounds = [i for i, x in enumerate(data) if "ATOM_TYPE" in x] + [len(data)]
    atom_blocks = [
        read_block_input(data[start:end], validator=defaultAPI)
        for start, end in zip(atom_type_bounds[:-1], atom_type_bounds[1:])
    ]

    return {
        "ion": {
            "atom_blocks": atom_blocks,
            "comments": new_comments,
            "sorting": {"sort": sort, "resort": resort},
        }
    }


@writer
def _write_aimd(
    fileobj,
    data_dict,
):
    raise NotImplementedError("Writing aimd file from python-api not supported!")
