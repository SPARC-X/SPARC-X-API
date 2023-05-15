#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

from warnings import warn


from setuptools.command.install import install
from warnings import warn


test_requires = [
    "pytest",
    "pyfakefs",
    "pytest-cov",
    "black",
    "flake8",
    "anybadge",
]

setup(
    name="sparc-dft-api",
    version="0.2",
    python_requires=">=3.8",
    description="Python Wrapper for the SPARC-X DFT Code",
    author="Ben Comer",
    author_email="ben.comer@gatech.edu",
    url="https://github.com/SPARC-X/sparc-dft-api",
    packages=find_packages(),
    install_requires=["ase>=3.22.0"],
    entry_points={
        "ase.io": [
            "sparc = sparc.io",
        ],
    },
    extras_require={
        "test": test_requires,
    },
    package_data={"sparc": ["psp/*", "sparc_json_api/*.json"]},
    include_package_data=True,
)
