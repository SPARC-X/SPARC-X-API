#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

from warnings import warn


from sparc.download_data import download_psp
try:
    download_psp()
except Exception:
    warn("Downloading external data to sparc-dft-api failed. Buidling wheel with them.")
    


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
    package_data={"sparc": ["psp/*"]},
    include_package_data=True,
)
