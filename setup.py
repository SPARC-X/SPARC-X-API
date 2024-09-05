#!/usr/bin/env python

from distutils.core import setup
from warnings import warn

from setuptools import find_packages
from setuptools.command.install import install

test_requires = [
    "pytest",
    "pyfakefs",
    "pytest-cov",
    "black",
    "isort",
    "flake8",
    "anybadge",
    "pre-commit",
]

setup(
    name="sparc-x-api",
    version="1.0.4",
    python_requires=">=3.8",
    description="Python API for the SPARC DFT Code",
    author="Tian Tian, Ben Comer",
    author_email="alchem0x2a@gmail.com, ben.comer@gatech.edu",
    url="https://github.com/SPARC-X/SPARC-X-API",
    packages=find_packages(),
    # ASE 3.22 dependency will be deprecated in 1.1.0+ release
    install_requires=["ase>=3.22.0", "numpy>=1.23", "packaging>=20.0", "psutil>=5.0.0"],
    entry_points={
        "ase.io": [
            "sparc = sparc.io",
        ],
        # The ioformats are only compatible with ase>=3.23
        "ase.ioformats": [
            "sparc = sparc.io:format_sparc",
            "sparc_ion = sparc.io:format_ion",
            "sparc_static = sparc.io:format_static",
            "sparc_geopt = sparc.io:format_geopt",
            "sparc_aimd = sparc.io:format_aimd",
        ],
        "console_scripts": ["sparc-ase=sparc.cli:main"],
    },
    extras_require={
        "test": test_requires,
    },
    package_data={"sparc": ["psp/*", "sparc_json_api/*.json"]},
    include_package_data=True,
)
