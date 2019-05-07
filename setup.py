#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='pysparc_x',
      version='0.1',
      description='Python Wrapper for the SPARC-X DFT Code',
      author='Ben Comer',
      author_email='ben.comer@gatech.edu',
      url='https://github.com/SPARC-X/pysparcx/tree/new_input',
      packages=find_packages(),
      install_requires=['spglib', 'numpy','ase','scipy'],
     )
