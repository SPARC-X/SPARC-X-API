#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='sparc-dft-api',
      version='0.1',
      python_requires='>3.5',
      description='Python Wrapper for the SPARC-X DFT Code',
      author='Ben Comer',
      author_email='ben.comer@gatech.edu',
      url='https://github.com/SPARC-X/sparc-dft-api',
      packages=find_packages(),
      install_requires=['spglib', 'numpy','ase','scipy'],
     )
