"""Provide a simple parser for pseudo potentials

The psp8 format is defined in abinit manual
https://docs.abinit.org/developers/psp8_info/

The first

"""

import os
import re
import shutil
from pathlib import Path
from warnings import warn

import numpy as np
from ase.data import atomic_names, chemical_symbols


class NotPSP8Format(Exception):
    def __init__(self, message):
        self.message = message


class NoMatchingPseudopotential(FileNotFoundError):
    def __init__(self, message):
        self.message = message


class MultiplePseudoPotentialFiles(Exception):
    def __init__(self, message):
        self.message = message


def parse_psp8_header(text):
    """Parse the first 4 lines of psp8 text
    if parsing failed, raise exception
    """
    header = text.split("\n")[:4]
    # Line 1

    psp8_magic = r"^\s*(?P<symbol>[\w]+)\s+ONCVPSP-(?P<psp8ver>[\d.\w]+)\s+r\_core=(?P<r_core>.*?)$"
    psp8_data = {}
    match = re.search(psp8_magic, header[0])
    if match is None:
        raise NotPSP8Format(f"The pseudopotential file is not in PSP8 format!")
    mgroup = match.groupdict()
    psp8_data["symbol"] = mgroup["symbol"].strip()
    psp8_data["psp8ver"] = mgroup["psp8ver"].strip()
    psp8_data["r_core"] = np.fromstring(mgroup["r_core"].strip(), sep=" ", dtype=float)
    # Line 2
    zatom, zion, pspd, *_ = header[1].split()
    psp8_data["zatom"] = float(zatom)
    psp8_data["zion"] = float(zion)
    # TODO: should we make date in datetime object?
    psp8_data["pspd"] = str(pspd).strip()

    # Line 3
    pspcod, pspxc, lmax, lloc, mmax, r2well, *_ = header[2].split()
    psp8_data["pspcod"] = int(pspcod)
    psp8_data["pspxc"] = int(pspxc)
    psp8_data["lmax"] = int(lmax)
    psp8_data["lloc"] = int(lloc)
    psp8_data["mmax"] = int(mmax)
    psp8_data["r2well"] = int(r2well)

    # Line 4
    rchrg, fchrg, qchrg, *_ = header[3].split()
    psp8_data["rchrg"] = float(rchrg)
    psp8_data["fchrg"] = float(fchrg)
    psp8_data["qchrg"] = float(qchrg)

    # Sanity check the symbol and zatom
    int_zatom = int(psp8_data["zatom"])
    if chemical_symbols[int_zatom] != psp8_data["symbol"]:
        raise NotPSP8Format(
            (
                f"The symbol defined in pseudo potential {psp8_data['symbol']} does not match "
                f"the Z={int_zatom}!"
            )
        )
    return psp8_data


def infer_pseudo_path(symbol, search_path):
    """Given an element symbol like 'Na', get the file name
    of the search_path (resolved) that search through the search path

    TODO: shall we support multiple directories?
    TODO: add a `setup` option like VASP?
    """
    search_path = Path(search_path).resolve()
    potfiles = (
        list(search_path.glob("*.psp8"))
        + list(search_path.glob("*.psp"))
        + list(search_path.glob("*.pot"))
    )
    candidates = []
    for pf in potfiles:
        try:
            psp8_data = parse_psp8_header(open(pf, "r").read())
        except Exception as e:
            print(e)
            psp8_data = None

        if psp8_data:
            if psp8_data["symbol"] == symbol:
                candidates.append(pf)
    if len(candidates) == 0:
        raise NoMatchingPseudopotential(
            (
                f"No pseudopotential file for {symbol} found "
                "under the search path {search_path}!"
            )
        )
    elif len(candidates) > 1:
        msg = (
            f"There are multiple psp8 files for {symbol}:\n"
            f"{candidates}. Please select the desired pseudopotential file!"
        )
        raise MultiplePseudoPotentialFiles(msg)

    else:
        return candidates[0]


def copy_psp_file(source_pot, target_dir, use_symbol=False):
    """Copy the pseudo potential file `source_pot` under `target_dir`

    if use_symbol is True, rename the potential to '{symbol}.psp8'

    the function returns the name of the pseudo potential file
    """

    source_pot = Path(source_pot)
    target_dir = Path(target_dir)
    psp8_data = parse_psp8_header(open(source_pot, "r").read())
    symbol = psp8_data["symbol"]
    if use_symbol:
        potname = f"{symbol}.psp8"
    else:
        potname = source_pot.name

    target_pot = target_dir / potname
    # shutil will copy
    shutil.copy(source_pot, target_pot)
    return potname


def find_pseudo_path(symbol, search_path=None, pseudopotential_mapping={}):
    """Get the pseudo potential file at best effort

    Searching priorities
    1) if pseudopotential_mapping has symbol as key, use the file name
       There are two possibilities
       i) filename does not contain directory information: i.e. Na-pbe.pot
          use search_path / filename for the mapping
       ii) filename contains directory information, directly use the file name
    2) No pseudopotential_mapping is given, get the psp from search_path
    """
    mapping_psp = pseudopotential_mapping.get(symbol, None)
    if mapping_psp is None:
        if search_path is None:
            raise NoMatchingPseudopotential(
                (
                    f"No psudopotentials found for {symbol} "
                    "because neither search_path nor psp name is provided."
                )
            )
        return infer_pseudo_path(symbol, search_path)
    else:
        str_psp = str(mapping_psp)
        mapping_psp = Path(mapping_psp)
        # if psp contains any path information (/, \\), treat is as a direct file
        is_node_file_name = (mapping_psp.name == str_psp) and (os.sep not in str_psp)
        if is_node_file_name:
            if search_path is None:
                raise NoMatchingPseudopotential(
                    (
                        f"You provide the pseudopotential name {mapping_psp} but no search path is defined. I cannot locate the pseudopotential file."
                    )
                )
            mapping_psp = Path(search_path) / str_psp

        if not mapping_psp.is_file():
            warn(
                (
                    f"Pseudopotential file {mapping_psp} is defined by user input but cannot be found!\n"
                    "Please check your setup. I'll write the .ion file anyway."
                )
            )
        return mapping_psp
