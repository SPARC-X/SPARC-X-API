#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 14:16:21 2018

@author: benjamin
"""

from ase import Atoms, Atom
from .utilities import valence_dict
from ase.units import Bohr
from ase.data import chemical_symbols, atomic_masses_iupac2016
import numpy as np
import os


def read_ion(fileobj):  # ,index):
    text = fileobj.read()
    comments_removed = []
    comments = []
    label = fileobj.name.split('.ion')[0]
        
    for line in text.split('\n'):
        comments_removed.append(line.split('#')[0])
        try:
            comments.append(line.split('#')[1])
        except:
            pass
    atoms = Atoms()
    ##################################
    # Parse the unit cell
    ##################################
    inpt_file_usable = True
    lat_vec_speficied = False
    comments_bad = False
    lat_array = []
    # Loop to read cell/latvec from either .inpt or the comment
    while True:
        # try to get the unit cell from the .inpt file in the same directory
        if label + '.inpt' in os.listdir('.') and inpt_file_usable == True:
            with open(label + '.inpt','r') as f:
                input_file = f.read()
            if 'CELL' not in input_file:
                inpt_file_usable = False
                del input_file
                continue 
            input_file = input_file.split('\n')
            for line in input_file:
                if 'CELL' in line:  # find the line with the cell info
                    cell = line.strip().split()[1:]
                    cell = np.array([float(a) for a in cell])
                if 'LATVEC' in line:  # get lattice vectors from next 3 lines
                    lat_vec_speficied = True
                    index = input_file.index(line)
                    for lat_vec in [input_file[a] for a in range(index + 1, index + 4)]:
                        vec = lat_vec.strip().split()
                        vec = np.array([float(a) for a in vec])
                        lat_array.append(vec / np.linalg.norm(vec))  # normalize
                    lat_array = np.array(lat_array)
            if lat_vec_speficied == False:
                lat_array = np.eye(3)
            if 'cell' in locals() and 'lat_array' in locals():
                atoms.cell = (lat_array.T * cell).T * Bohr
                break  # we got the cell, leave the while loop
            else:
                inpt_file_usable = False
                del input_file
                continue
        
        elif comments != [] and comments_bad == False:
            if 'CELL' in comments[1]:  # Check only the second line for a CELL
                cell = np.empty((3,3))
                try:
                    cell = comments[1].strip().split()[1:]
                    cell = [float(a) for a in cell]
                    for lat_vec in comments[3:6]:  # check only these lines
                        vec = lat_vec.strip().split()
                        vec = np.array([float(a) for a in vec])
                        lat_array.append(vec / np.linalg.norm(vec))  # normalize
                    lat_array = np.array(lat_array)
                    atoms.cell = (lat_array.T * cell).T * Bohr
                    break
                except:
                    comments_bad = True
            else:  # if getting it from the comments fails, return 0 unit cell
                raise(Warning('No lattice vectors were found in either the .inpt file or in the comments of the .ion file. Thus no unit cell was set for the resulting atoms object. Use the output atoms at your own risk'))
                atoms.cell = np.zeros((3, 3))
                break    
        else:  # if there is no cell in the .inpt file, and the .ion file, return 0 unit cell
            raise(Warning('No lattice vectors were found in either the .inpt file or in the comments of the .ion file. Thus no unit cell was set for the resulting atoms object. Use the output atoms at your own risk'))
            atoms.cell = np.zeros((3,3))
            break
    ############################################
    # parse the atoms
    ############################################
    for i, line in enumerate(comments_removed):
        if 'ATOM_TYPE:' in line:
            element = line.strip().split()[1]
            num = int(comments_removed[i + 1].strip().split()[1])
            pseudo = comments_removed[i + 2].strip().split()[1]
            frac = 'COORD_FRAC' in comments_removed[i + 4]
            if frac and sum(sum(atoms.cell)) == 0:
                raise(Exception('There is no specified unit cell and the coordinates are in fractional units. This file cannot be converted to an atoms object'))
            if 'ATOMIC MASS' in ''.join(comments_removed) or \
               'ATOMIC_MASS' in ''.join(comments_removed):
                offset = 5
            else:
                offset = 4
            for coord_set in comments_removed[i + offset: i + num + offset]:
                if frac:
                    x1, x2, x3 = [float(a) for a in coord_set.split()]
                    x, y, z = sum([x * a for x, a in zip([x1, x2, x3], atoms.cell)])
                else:
                    x, y, z = [float(a) * Bohr for a in coord_set.split()]
                atoms += Atom(symbol = element, position=(x, y, z))
    return atoms


def write_ion(fileobj, atoms, pseudo_dir = None, scaled = True, comment = ''):
    """
    Standard ase io file for reading the sparc-x .ion format
    
    inputs:
        atoms (ase atoms object):
            an ase atoms object of the system being written to a file
        pseudos (list):
            a list of the locations of the pseudopotential files to be used
    """
    elements = sorted(list(set(atoms.get_chemical_symbols())))
    
    fileobj.write('# Input File Generated By SPARC ASE Calculator #\n')
    fileobj.write('#CELL:')
    for comp in atoms.cell:
        fileobj.write(' ' + str(np.linalg.norm(comp) /  Bohr))
    fileobj.write('\n#LATVEC\n')
    for comp in atoms.cell:
        fileobj.write('#')
        comp_n = comp / np.linalg.norm(comp)  # normalize
        for i in comp_n:
            fileobj.write(' ' + str(i))
        fileobj.write('\n')
    fileobj.write('# ' + comment + '\n\n\n')
    for i, element in enumerate(elements):
        fileobj.write('ATOM_TYPE: ')
        fileobj.write(element + '\n')

        fileobj.write('N_TYPE_ATOM: ')
        fileobj.write(str(atoms.get_chemical_symbols().count(element)) + '\n')

        if pseudo_dir is not None:
            pseudos_in_dir = [a for a in os.listdir(os.environ['SPARC_PSP_PATH']) \
                        if a.endswith(element+'.pot')]
            filename = [a for a in os.listdir(os.environ['SPARC_PSP_PATH']) \
                        if a.endswith(element+'.pot')][0]
            os.system('cp $SPARC_PSP_PATH/' + filename + ' .')
            fileobj.write('PSEUDO_POT: {}\n'.format(filename)) 
            
            atomic_number = chemical_symbols.index(element)
            atomic_mass = atomic_masses_iupac2016[atomic_number]
            fileobj.write('ATOMIC_MASS: {}\n'.format(atomic_mass))
        
        if scaled == False:
            fileobj.write('COORD:\n')
            positions = atoms.get_positions()
            positions /= Bohr
        else:
            fileobj.write('COORD_FRAC:\n')
            positions = atoms.get_scaled_positions()
        for atom, position in zip(atoms, positions):
            if atom.symbol == element:
                for component in position:
                    fileobj.write('    ' + str(component))
                fileobj.write('\n')
        fileobj.write('\n\n')
