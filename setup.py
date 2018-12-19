import sys
import os
script_path = os.path.dirname(os.path.abspath( __file__ ))



print('This is not a sophisticated installation file. To install this python package you should add the\ndirectory this script is in to your PYTHONPATH environment varible. here\'s a suggested line to add\nto your ~/.bashrc file:')
print('export PYTHONPATH=' + script_path + ':$PYTHONPATH')

"""
from setuptools import setup

setup(name='pysparcx',
          version='0.1',
          description='The python wrapper for the ab inito code SPARC-X',
          url='https://github.com/SPARC-X/ase-sparc.git',
          author='Flying Circus',
          author_email='',
          license='MIT',
          install_requires=[
          'ase',
          'numpy',
          'scipy'],
          packages=['pysparcx'],
          zip_safe=False)
"""
