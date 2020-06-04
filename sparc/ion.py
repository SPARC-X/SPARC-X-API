# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 14:16:21 2018

Ben Comer (Georgia Tech)
"""
import shutil
import os

import numpy as np

from ase import Atoms, Atom
from ase.units import Bohr
from ase.data import chemical_symbols, atomic_masses_iupac2016
import warnings


def read_ion(fileobj, recover_indices=True, recover_constraints=True):
    text = fileobj.read()
    comments_removed = []
    comments = []
    label = fileobj.name.split('.ion')[0]

    for line in text.splitlines():  # break into lines
        # remove and store the comments
        entry = line.split('#')
        if not entry[0]:
            pass
        else:
            comments_removed.append(entry[0].strip())
        try:
            comments.append(line.split('#')[1])
        except:
            pass
    atoms = Atoms()
    ##################################
    # Parse the unit cell
    ##################################

    """
    because the unit cell is not included in the .ion atomic positions
    file in SPARC, this interface writes the information into the comments
    of .ion files. If a .inpt file is present this code will read that
    if it is not it will attempt to read the comments to find that information.
    The logic for doing this is quite complicated (as seen below)
    """

    inpt_file_usable = True
    lat_vec_speficied = False
    comments_bad = False
    lat_array = []
    # Loop to read cell/latvec from either .inpt or the comment
    # this is pretty complicated
    while True:
        # try to get the unit cell from the .inpt file in the same directory
        if label + '.inpt' in os.listdir('.') and inpt_file_usable == True:
            with open(label + '.inpt', 'r') as f:
                input_file = f.read()
            if 'CELL' not in input_file:
                # We can't find the CELL in the input file
                # set the flag to false and re-run the while loop.
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
                        lat_array.append(
                            vec / np.linalg.norm(vec))  # normalize
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

        # if the input file isn't usable, check the comments of the .ion file
        elif comments != [] and comments_bad == False:
            if 'CELL' in comments[1]:  # Check only the second line for a CELL
                cell = np.empty((3, 3))
                try:
                    cell = comments[1].strip().split()[1:]
                    cell = [float(a) for a in cell]
                    for lat_vec in comments[3:6]:  # check only these lines
                        vec = lat_vec.strip().split()
                        vec = np.array([float(a) for a in vec])
                        lat_array.append(
                            vec / np.linalg.norm(vec))  # normalize
                    lat_array = np.array(lat_array)
                    atoms.cell = (lat_array.T * cell).T * Bohr
                    break
                except:
                    comments_bad = True
            else:  # if getting it from the comments fails, return 0 unit cell
                warnings.warn('No lattice vectors were found in either the .inpt'
                              ' file or in the comments of the .ion file. Thus no'
                              ' unit cell was set for the resulting atoms object. '
                              'Use the output atoms at your own risk')
                atoms.cell = np.zeros((3, 3))
                break
        else:  # if there is no cell in the .inpt file, and the .ion file, return 0 unit cell
            warnings.warn('No lattice vectors were found in either the .inpt file '
                          'or in the comments of the .ion file. Thus no unit cell '
                          'was set for the resulting atoms object. Use the output '
                          'atoms at your own risk')
            atoms.cell = np.zeros((3, 3))
            break

    ############################################
    # parse the atoms
    ############################################

    # parse apart the comments to try to recover the indices and boundary conditions

    """
    The strategy of this code is to get the input text separated from the
    comments.

    The comments are then used to glean the information that is
    normally stored in the .inpt file, as well as the original indices
    of the atoms if this file was made by this wrapper.

    from there figure out where the different "Atom Type" blocks are 
    located in the full text with the comments removed. Once that has
    been found parse these sections to gain recover the atomic positions
    and elemental identies of these "atom types."

    We also need to find the locations of the "RELAX" blocks that contain
    the information on which atoms are constained in each "Atom Type"
    block.
    """

    indices_from_comments = []
    constraints = []
    spins = []
    for comment in comments:
        if 'index' in comment:
            if len(comment.split()) == 2:
                try:
                    index = int(comment.split()[1])
                    indices_from_comments.append(index)
                except:
                    pass
        if 'PBC:' in comment:
            pbc_list = []
            pbc = comment.split()[1:]
            for c in pbc:
                if c == 'True' or c == 'true':
                    pbc_list.append(True)
                else:
                    pbc_list.append(False)
            atoms.set_pbc(pbc_list)
            del pbc_list, pbc
    # find the index of line for all the different atom types
    atom_types = [i for i, x in enumerate(
        comments_removed) if 'ATOM_TYPE:' in x]
    relax_blocks = [i for i, x in enumerate(comments_removed) if 'RELAX:' in x]
    spin_blocks = [i for i, x in enumerate(comments_removed) if 'SPIN:' in x]
    for i, atom_type in enumerate(atom_types):
        type_dict = {}
        if i == len(atom_types) - 1:  # treat the last block differently
            # Get the slice of text associated with this atom type
            type_slice = comments_removed[atom_types[i]:]

            # figure out if there are constraints after this block
            if recover_constraints:
                relax_block_index = [a for a in relax_blocks if a > atom_type]
            spin_block_index = [a for a in spin_blocks if a > atom_type]
        else:
            # Get the slice of text associated with this atom type
            type_slice = comments_removed[atom_types[i]: atom_types[i+1]]
            [a for a in relax_blocks if a > atom_type]

            # figure out if there are constraints after this block.
            # the constraint index will be sandwiched between the indicies
            # the current block and the next block
            if recover_constraints:
                relax_block_index = [
                    a for a in relax_blocks if a > atom_types[i]]
                relax_block_index = [
                    a for a in relax_block_index if a < atom_types[i+1]]
            spin_block_index = [a for a in spin_blocks if a > atom_types[i]]
            spin_block_index = [
                a for a in spin_block_index if a < atom_types[i+1]]

        # extract informaton about the atom type from the section header
        for info in ['PSEUDO_POT', 'ATOM_TYPE', 'ATOMIC_MASS', 'COORD',
                     'N_TYPE_ATOM']:
            for line in type_slice[:15]:  # narrow the search for speed
                if info in line:
                    if 'COORD' in line:
                        if 'FRAC' in line:
                            type_dict['COORD_FRAC'] = 1
                        else:
                            type_dict['COORD'] = 1
                    elif 'COORD' not in line:
                        type_dict[info] = line.split()[1]

        # get the lines that contain the constraints block
        if recover_constraints:
            if len(relax_block_index) == 0:
                pass
            elif len(relax_block_index) == 1:
                # offest by one line
                relax_block_index = relax_block_index[0] + 1
                relax_block_end = relax_block_index + \
                    int(type_dict['N_TYPE_ATOM'])
                relax_slice = comments_removed[relax_block_index: relax_block_end]
            elif len(relax_block_index) > 1:
                raise Exception('There appear to be multiple blocks of'
                                ' constraints in one or more of the atom'
                                ' types in your .ion file. Please inspect'
                                ' it to repair it or pass in '
                                '`recover_constraints = False` to ingore'
                                ' constraints')
        # the same as the code above, but for spin
        if len(spin_block_index) == 0:
            pass
        elif len(spin_block_index) == 1:
            spin_block_index = spin_block_index[0] + 1  # offest by one line
            spin_block_end = spin_block_index + int(type_dict['N_TYPE_ATOM'])
            spin_slice = comments_removed[spin_block_index: spin_block_end]
        elif len(spin_block_index) > 1:
            raise Exception('There appear to be multiple blocks of'
                            ' spin values in one or more of the atom'
                            ' types in your .ion file. Please inspect'
                            ' it to repair it or pass in')

        # now parse out the atomic positions
        for coord_set in type_slice[len(type_dict):int(type_dict['N_TYPE_ATOM']) + len(type_dict)]:
            if 'COORD_FRAC' in type_dict.keys():
                x1, x2, x3 = [float(a) for a in coord_set.split()[:3]]
                x, y, z = sum(
                    [x * a for x, a in zip([x1, x2, x3], atoms.cell)])
            elif 'COORD' in type_dict.keys():
                x, y, z = [float(a) * Bohr for a in coord_set.split()[:3]]
            atoms += Atom(symbol=type_dict['ATOM_TYPE'], position=(x, y, z))
        # get the constraints
        if recover_constraints:
            if 'relax_slice' in locals():
                for cons_set in relax_slice:
                    constraints.append([int(a) for a in cons_set.split()])
                del relax_slice
            else:
                # there aren't constraints with this block, put in empty lists
                constraints += [[]] * int(type_dict['N_TYPE_ATOM'])
        if 'spin_slice' in locals():
            for init_spin in spin_slice:
                spins.append(float(init_spin))
            del spin_slice
        else:
            # there aren't spins with this block, put in zeros
            spins += [0] * int(type_dict['N_TYPE_ATOM'])
    # check if we can reorganize the indices
    if len(indices_from_comments) == len(atoms) and recover_indices:
        new_atoms = Atoms(['X'] * len(atoms),
                          positions=[(0, 0, 0)] * len(atoms))
        new_atoms.set_cell(atoms.cell)
        new_spins = [None] * len(atoms)
        # reassign indicies
        for old_index, new_index in enumerate(indices_from_comments):
            new_atoms[new_index].symbol = atoms[old_index].symbol
            new_atoms[new_index].position = atoms[old_index].position
            new_atoms.pbc = atoms.pbc
            new_spins[new_index] = spins[old_index]
        assert new_atoms.get_chemical_formula() == atoms.get_chemical_formula()
        atoms = new_atoms
        spins = new_spins
        # reorganize the constraints now
        if recover_constraints:
            new_constraints = [0] * len(atoms)
            for old_index, new_index in enumerate(indices_from_comments):
                new_constraints[new_index] = constraints[old_index]
            constraints = new_constraints

    atoms.set_initial_magnetic_moments(spins)
    # add constraints
    if recover_constraints:
        constraints = decipher_constraints(constraints)
        atoms.set_constraint(constraints)
    return atoms


def decipher_constraints(constraints):
    """
    This is just a helper function to translate from a list of
    constraints from SPARC to a set of ASE constraints

    Parameters:
        constraints (list):
            a list for lists, either empty or in the form [x,y,z]
            i.e. [[],[0,0,1],[1,1,1]]

    returns:
        cons_list (list)
            a list of ase contraint classes
    """
    import ase.constraints as c
    cons_list = []
    fix_atoms = []
    for i, constraint in enumerate(constraints):
        if constraint.count(1) == 3 or constraint == []:
            continue
        elif constraint.count(1) == 0:
            fix_atoms.append(i)
        elif constraint.count(1) == 1:
            cons_list.append(c.FixedLine(i, constraint))
        elif constraint.count(1) == 2:
            tmp = [int(not a) for a in constraint]  # flip ones and zeros
            cons_list.append(c.FixedPlane(i, tmp))
    cons_list.append(c.FixAtoms(fix_atoms))
    return cons_list


def write_ion(fileobj, atoms, pseudo_dir=None, scaled=True,
              add_constraints=True, copy_psp=True, comment=''):
    """
    Standard ase io file for reading the sparc-x .ion format

    inputs:
        atoms (ase atoms object):
            an ase atoms object of the system being written to a file
        pseudos (list):
            a list of the locations of the pseudopotential files to be used
    """

    directory = os.path.dirname(fileobj.name)

    elements = sorted(list(set(atoms.get_chemical_symbols())))

    fileobj.write('# Input File Generated By SPARC ASE Calculator #\n')
    fileobj.write('#CELL:')
    for comp in atoms.cell:
        fileobj.write('  {}'.format(
            format(np.linalg.norm(comp) / Bohr, ' .15f')))
    fileobj.write('\n#LATVEC\n')
    for comp in atoms.cell:
        fileobj.write('#')
        comp_n = comp / np.linalg.norm(comp)  # normalize
        for i in comp_n:
            fileobj.write(' {}'.format(format(i, ' .15f')))
        fileobj.write('\n')
    if atoms.pbc is not None:
        fileobj.write('#PBC: ')
        for condition in atoms.pbc:
            fileobj.write(str(condition) + ' ')
        fileobj.write('\n')
    fileobj.write('# ' + comment + '\n\n\n')

    # translate the constraints to usable strings
    if add_constraints == True:
        cons_indices = []
        cons_strings = []
        for constraint in atoms.constraints:
            name = type(constraint).__name__
            if name == 'FixAtoms':
                cons_indices += list(constraint.index)
                cons_strings += ['0 0 0\n'] * len(constraint.index)
            elif name == 'FixedLine' or 'FixedPlane':
                line_dir = [int(np.ceil(a)) for a in constraint.dir]
                max_dirs = 1  # for FixedLine
                # The API of FixedLine and FixedPlane is the same
                # except for the need for an inversion of indices
                if name == 'FixedPlane':
                    line_dir = [int(not a) for a in line_dir]
                    max_dirs = 2
                if line_dir.count(1) > max_dirs:
                    warnings.warn('The .ion filetype can only support'
                                  'freezing entire dimensions (x,y,z).'
                                  'The {} constraint on this atoms'
                                  ' Object is being converted to allow '
                                  'movement in any directions in which the'
                                  ' atom had some freedom to move'.format(name))
                cons_indices += constraint.get_indices()
                cons_strings.append('{} {} {}\n'.format(*line_dir))
            else:
                warnings.warn('The constraint type {} is not supported by'
                              ' the .ion format. This constraint will be'
                              ' ignored')

        # just make sure that there aren't repeats
        assert len(set(cons_indices)) == len(cons_indices)

    for i, element in enumerate(elements):
        fileobj.write('ATOM_TYPE: ')
        fileobj.write(element + '\n')

        fileobj.write('N_TYPE_ATOM: ')
        fileobj.write(str(atoms.get_chemical_symbols().count(element)) + '\n')
        if add_constraints:
            constraints_string = 'RELAX:\n'
        magmoms = [float(a) for a in atoms.get_initial_magnetic_moments()]
        if magmoms != [0.] * len(atoms):
            spin_string = 'SPIN:\n'

        # TODO: fix this psuedopotential finding code
        if pseudo_dir is not None:
            if not os.path.isdir(pseudo_dir):
                pseudo_dir = os.environ['SPARC_PSP_PATH']
                warnings.warn('The path entered for `pseudo_dir` could '
                              'not be found. Using SPARC_PSP_PATH and generic pseudopotential '
                              'file names (<symbol>.pot). Psuedopotentials must be added '
                              'manually for SPARC to execute successfully.')
                fileobj.write('PSEUDO_POT: {}.pot\n'.format(element))
            else:
                pseudos_in_dir = [a for a in os.listdir(pseudo_dir)
                                  if a.endswith(element+'.pot')]

                filename = [a for a in os.listdir(pseudo_dir)
                            if a.endswith(element+'.pot')]

                if len(filename) == 0:
                    filename = element+'.pot'
                    warnings.warn('No pseudopotential detected for '+element+'. '
                                  'Using generic filename ('+element +
                                  '.pot). The pseudopotential '
                                  ' must be added manually for SPARC to execute successfully.')
                else:
                    if len(filename) > 1:
                        warnings.warn('Multiple psudopotentials detected for '
                                      '{} ({}).'.format(element, str(filename))+' Using '+filename[0]+'.')
                    if copy_psp:
                        filename = filename[0]
                        full_psp_path = os.path.abspath(os.path.join(pseudo_dir, filename))
                        full_dest_path = os.path.abspath(os.path.join(directory, filename))
                        if full_psp_path != full_dest_path:
                            #shutil.copyfile(os.path.join(pseudo_dir, filename),
                            #                os.path.join(directory, filename))
                            shutil.copyfile(full_psp_path, full_dest_path)
                    else:
                        filename = os.path.join(pseudo_dir, filename[0])

                #os.system('cp $SPARC_PSP_PATH/' + filename + ' .')
                fileobj.write('PSEUDO_POT: {}\n'.format(filename))

        else:
            fileobj.write('PSEUDO_POT: {}.pot\n'.format(element))

        atomic_number = chemical_symbols.index(element)
        atomic_mass = atomic_masses_iupac2016[atomic_number]
        fileobj.write('ATOMIC_MASS: {}\n'.format(atomic_mass))

        if scaled == False:
            fileobj.write('COORD:\n')
            positions = atoms.get_positions(wrap=False)
            positions /= Bohr
        else:
            fileobj.write('COORD_FRAC:\n')
            positions = atoms.get_scaled_positions(wrap=False)
        for atom, position in zip(atoms, positions):
            if atom.symbol == element:
                for component in position:
                    fileobj.write('    {}'.format(format(component, ' .15f')))
                fileobj.write('   # index {}'.format(atom.index))
                fileobj.write('\n')
                # mess with the constraints
                if add_constraints:
                    if atom.index in cons_indices:
                        constraints_indices_index = cons_indices.index(
                            atom.index)
                        constraints_string += cons_strings[constraints_indices_index]
                    else:
                        constraints_string += '1 1 1\n'
                if 'spin_string' in locals():
                    spin_string += format(atom.magmom, ' .15f') + '\n'
        # dump in the constraints
        if add_constraints:
            fileobj.write(constraints_string)
        if 'spin_string' in locals():
            fileobj.write(spin_string)
        fileobj.write('\n\n')
