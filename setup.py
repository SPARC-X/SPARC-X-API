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

doc_requires = [
    "sphinx_rtd_theme",
    "sphinx_tabs",
    "myst-parser",
]

setup(
    name="sparc-x-api",
    version="1.0.5",
    python_requires=">=3.8",
    description="Python API for the SPARC DFT Code",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Tian Tian, Lucas R Timmerman, Ben Comer",
    author_email="alchem0x2a@gmail.com, ltimmerman3@gatech.edu, ben.comer@gatech.edu",
    url="https://github.com/SPARC-X/SPARC-X-API",
    packages=find_packages(),
    # ASE 3.22 dependency will be deprecated in 2.0+ release
    install_requires=["ase>=3.22.0", "numpy>=1.23", "packaging>=20.0", "psutil>=5.0.0"],
    entry_points={
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
        "doc": test_requires + doc_requires,
    },
    package_data={"sparc": ["psp/*.psp8", "sparc_json_api/*.json"]},
    include_package_data=True,
)
