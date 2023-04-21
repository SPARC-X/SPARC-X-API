#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

import os
from pathlib import Path
import tempfile
import urllib.request
import zipfile
from setuptools import Command, find_packages, setup
from io import BytesIO
import shutil

repo_root = Path(__file__).parent

class DownloadPSPs(Command):
    """Download the SPSM pseudopotentials associated with SPARC"""

    description = "Download and extract pseudopotentials distributed by SPARC"
    user_options = []
    download_url = "https://github.com/SPARC-X/SPARC/archive/b702c1061400a2d23c0e223e32182609d7958156.zip"

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Download the ZIP archive of public SPARC (last release of PSPs)
        url = self.download_url
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with urlopen(download_url) as zipresp:
                with zipfile.ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall(tmpdir)
                    source_dir = tmpdir / "SPARC-b702c1061400a2d23c0e223e32182609d7958156" / "psp"
                    if not source_dir.is_dir():
                        raise FileNotFoundError("Error downloading or extracting zip")
                    shutil.copy(source_dir, repo_root / "psp")


test_requires = ["pytest", "pyfakefs", "pytest-cov", "black", "flake8", "anybadge"]

setup(
    name="sparc-dft-api",
    version="0.2",
    python_requires=">=3.8",
    description="Python Wrapper for the SPARC-X DFT Code",
    author="Ben Comer",
    author_email="ben.comer@gatech.edu",
    url="https://github.com/SPARC-X/sparc-dft-api",
    packages=find_packages(),
    install_requires=["spglib", "numpy>=1.20", "ase>=3.22.0", "scipy"],
    entry_points={
        "ase.io": [
            "sparc = sparc.sparc_io_bundle",
        ],
    },
    extras_require={
        "test": test_requires,
    },
    # cmdclass={
    #     "download_sparc_psps": DownloadPSPs,
    # },
    package_data={"sparc": ["../psp/"]},
)
