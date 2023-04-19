"""Provide a simple parser for pseudo potentials

The psp8 format is defined in abinit manual
https://docs.abinit.org/developers/psp8_info/

The first 

"""

import re
import os
import numpy as np
from pathlib import Path
from ase.data import atomic_names, chemical_symbols


class NotPSP8Format(Exception):
    def __init__(self, message):
        self.message = message


def parse_psp8_header(text):
    """Parse the first 4 lines of psp8 text
    if parsing failed, raise exception
    """
    header = text.split("\n")[:4]
    # Line 1
    
    psp8_magic = (
    r"^\s*(?P<symbol>[\w]+)\s+ONCVPSP-(?P<psp8ver>[\d.\w]+)\s+r\_core=(?P<r_core>.*?)$")
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
        raise NotPSP8Format((
            f"The symbol defined in pseudo potential {psp8_data['symbol']} does not match "
            f"the Z={int_zatom}!"
        ))
    return psp8_data

    
    
