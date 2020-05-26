#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 12:52:34 2018

@author: benjamin
"""

import numpy as np
from ase.data import chemical_symbols
from ase.units import Bohr, Hartree, fs, GPa
from ase.io.jsonio import encode
from ase.atoms import Atoms
from ase.atom import Atom
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.trajectory import Trajectory
#from scipy.misc import factorial
from collections import OrderedDict
import json
import os

valences = [0]+[1, 2]+[1, 2, 3, 4, 5, 6, 7, 8]*2 + \
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]*2 + \
    list(np.arange(1, 33))*2
# atomic numbers correspond to traditional number of valence electrons
valence_dict = {}
for number, element in enumerate(chemical_symbols):
    valence_dict[element] = valences[number]

"""
def h2gpts(h, cell_cv, idiv=4):
"""
"""Convert grid spacing to number of grid points divisible by idiv.
    Taken from GPAW:
        https://gitlab.com/gpaw/gpaw/blob/master/gpaw/utilities/__init__.py

    Note that units of h and cell_cv must match!

    h: float
        Desired grid spacing in.
    cell_cv: 3x3 ndarray
        Unit cell.
"""
"""
    L_c = (np.linalg.inv(cell_cv)**2).sum(0)**-0.5
    return np.maximum(idiv, (L_c / h / idiv + 0.5).astype(int) * idiv)
"""


def h2gpts(h, cell_cv, idiv=4):
    cell_lengths = np.linalg.norm(cell_cv, axis=1)
    grid = np.ceil(cell_lengths/h)
    grid = np.maximum(idiv, grid)
    return [int(a) for a in grid]


def atoms_dict(atoms):
    d = OrderedDict(atoms=[{'symbol': atom.symbol,
                            # 'position': json.loads(encode(atom.position)),
                            'position': [float(a) for a in atom.position],
                            'tag': atom.tag,
                            'index': atom.index,
                            'charge': atom.charge,
                            # 'momentum': json.loads(encode(atom.momentum)),
                            'momentum': [float(a) for a in atom.momentum],
                            'magmom': atom.magmom}
                           for atom in atoms],
                    cell=list(atoms.cell.array),
                    pbc=atoms.pbc,
                    info=atoms.info,
                    constraints=[c.todict() for c in atoms.constraints])
    return d


def dict_atoms(d):
    atoms = Atoms([Atom(atom['symbol'],
                        position=atom['position'],
                        tag=atom['tag'],
                        momentum=atom['momentum'],
                        magmom=atom['magmom'],
                        charge=atom['charge'])
                   for atom in d['atoms']],
                  cell=d['cell'],
                  pbc=d['pbc'],
                  info=d['info'],
                  constraint=[dict2constraint(c) for c in d['constraints']])
    return atoms


def parse_output(label='sprc-calc', calc_type=None, write_traj=False):
    """
    Parses almost all useful information from the SPARC
    output file

    inputs:
        label (str):
            The base name of the output files

        write_traj (bool):
            If set to True, a trajectory file will be written
            from the outputfile.

    returns:
        atoms (ASE atoms object):
            An ASE atoms object from the last step. This is mainly
            used to set the new atoms for the calculator after 
            running SPARC using the internal relaxation algorithms.

        input_dict (dict):
            A dictionary containing the original input arguements
            put into SPARC. This can be used to rebuild a calculator
            object (e.g. SPARC(atoms = atoms, **input_dict))
    """
    f = open(label + '.out', 'r')
    text = f.read()
    f.close()
    text = text.split('SPARC')[-1]  # Get only the last run in the file

    # Recover the original input parameters
    input_parameters = text.split('*' * 75)[2]  # Use this later
    input_parameters = input_parameters.split('\n')
    input_dict = {}
    in_lattice_block = False
    for input_arg in input_parameters[1:-1]:
        # once we find LATVEC, we analyze the next 3 lines differently
        if 'LATVEC' in input_arg or in_lattice_block:
            if not 'blk_ln' in locals():
                input_dict['LATVEC'] = []
                in_lattice_block = True
                blk_ln = 1
                continue
            lat_vec = [float(a) for a in input_arg.strip().split()]
            input_dict['LATVEC'].append(lat_vec)
            if blk_ln == 3:
                in_lattice_block = False
                input_dict['LATVEC'] = np.array(input_dict['LATVEC'])
                continue
            else:
                blk_ln += 1
                continue

        # print(input_arg)
        kw, arg = input_arg.split(':')
        input_dict[kw.strip()] = arg.strip()
        if len(arg.split()) > 1:  # Some arugments are lists
            input_dict[kw.strip()] = arg.split()
    input_dict['label'] = input_dict['OUTPUT_FILE']
    del input_dict['OUTPUT_FILE']

    cell = [float(a) for a in input_dict['CELL']]
    cell = np.eye(3) * cell * Bohr
    if input_dict['BOUNDARY_CONDITION'] == '2':
        pbc = [True, True, True]
    else:
        pbc = [False, False, False]

    # Figure out how many 'types' of atoms there are
    s = os.popen('grep "Total number of atom types" ' + label + '.out')
    ntypes = int(s.readlines()[-1].split(':')[1].strip())
    s.close()

    # For each type, grep out the chemical symbol
    s = os.popen('grep "Atom type" ' + label + '.out')
    atom_elements = s.readlines()[-ntypes:]
    elements = [a.split()[-2] for a in atom_elements]
    s.close()

    # For each type, grep out how many of that type there are
    s = os.popen('grep "Number of atoms of type" ' + label + '.out')
    num_elements = s.readlines()[-ntypes:]
    numbers = [int(a.split()[-1]) for a in num_elements]
    s.close()

    # Make a list containing the elements of each atom in order
    chemical_symbols = []
    for sym, num in zip(elements, numbers):
        chemical_symbols += num * [sym]

    # if not os.path.isfile(label + '.relax'):
    # if os.path.isfile(label + '.relax'):
    if calc_type == 'relax':
        return parse_relax(label, write_traj=write_traj,
                           pbc=pbc, cell=cell,
                           chemical_symbols=chemical_symbols), input_dict
    elif calc_type == 'MD':
        return parse_MD(label, write_traj=write_traj,
                        pbc=pbc, cell=cell,
                        chemical_symbols=chemical_symbols), input_dict


def parse_relax(label, write_traj=False,
                pbc=False, cell=None, chemical_symbols=[]):
    f = open(label + '.relax')
    #f = open(label + '.restart')
    text = f.read()
    # Parse out the steps
    if text == '':
        return None
    steps = text.split(':RELAXSTEP:')[1:]
    if write_traj == False:
        steps = [steps[-1]]
    else:
        traj = Trajectory(label + '.traj', mode='w')

    # Parse out the energies
    n_geometric = len(steps)
    s = os.popen('grep "Total free energy" ' + label + '.out')
    engs = s.readlines()[-n_geometric:]
    engs = [float(a.split()[-2]) * Hartree for a in engs]
    s.close()

    # build a traj file out of the steps
    for j, step in enumerate(steps):
        positions = step.split(':')[2].strip().split('\n')
        forces = step.split(':')[4].strip().split('\n')
        frc = np.empty((len(forces), 3))
        atoms = Atoms()
        for i, f in enumerate(forces):
            frc[i, :] = [float(a) * Hartree / Bohr for a in f.split()]
            atoms += Atom(chemical_symbols[i],
                          [float(a) * Bohr for a in positions[i].split()])

        atoms.set_calculator(SinglePointCalculator(atoms, energy=engs[j],
                                                   forces=frc))
        atoms.set_pbc(pbc)
        atoms.cell = cell
        if write_traj == True:
            traj.write(atoms)
    atoms.set_calculator()
    return atoms


def parse_MD(label, write_traj=False, pbc=False, cell=None, chemical_symbols=[]):
    f = open(label + '.aimd')
    #f = open(label + '.restart')
    text = f.read()
    if text == '':
        return None
    # Parse out the steps
    steps = text.split(':MDSTEP:')[1:]
    if write_traj == False:
        steps = [steps[-1]]
    else:
        traj = Trajectory(label + '.traj', mode='w')

    # Parse out the energies
    n_images = len(steps)
    s = os.popen('grep ":FEN:" ' + label + '.aimd')
    engs = s.readlines()[-n_images:]
    engs = [float(a.split()[-1]) * Hartree for a in engs]
    s.close()

    # build a traj file out of the steps
    for j, step in enumerate(steps):
        # Find Indicies
        colons = step.split(':')
        #pos_index = colons.index('R(Bohr)') + 1
        #frc_index = colons.index('F(Ha/Bohr)') + 1
        #vel_index = colons.index('V(Bohr/atu)') + 1
        pos_index = colons.index('R') + 1
        frc_index = colons.index('F') + 1
        vel_index = colons.index('V') + 1
        # Parse the text
        positions = colons[pos_index].strip().split('\n')
        forces = colons[frc_index].strip().split('\n')
        velocities = colons[vel_index].strip().split('\n')
        # Initialize the arrays
        frc = np.empty((len(forces), 3))
        vel = np.empty((len(velocities), 3))
        stress = np.zeros((3, 3))
        atoms = Atoms()
        for i, f, v in zip(range(len(forces)), forces, velocities):
            frc[i, :] = [float(a) * Hartree / Bohr for a in f.split()]
            vel[i, :] = [float(a) / Bohr / fs for a in v.split()]
            atoms += Atom(chemical_symbols[i],
                          [float(a) * Bohr for a in positions[i].split()])
        if 'STRESS' in step:
            stress_index = colons.index('STRESS_TOT(GPa)') + 1
            for i, s in enumerate(colons[stress_index].strip().split('\n')):
                stress[i, :] = [float(a) * GPa for a in s.split()]
        atoms.set_velocities(vel)
        atoms.set_pbc(pbc)
        atoms.cell = cell
        atoms.set_calculator(SinglePointCalculator(atoms,
                                                   energy=engs[j] * len(atoms),
                                                   stress=stress,
                                                   forces=frc))
        if write_traj == True:
            traj.write(atoms)
    atoms.set_calculator()
    return atoms


def Ecut2h(FDn, Ecut, tol=0.1):
    # to do complete python conversion and integrate
    """
    converted from MATLAB code provided by Qimen
    """
    w2 = np.zeros(FDn+1)
    for i in range(FDn+1):
        k = i+1
        w2[i+1] = (2*(-1)**(k+1))*(factorial(FDn)**2) / \
            (k*k*factorial(FDn-k)*factorial(FDn+k))
        w2[0] = w2[0]-2*(1/(k*k))
    kk = np.linspace(0, np.pi, 1001)
    y_cos = -w2[0] + (-2*w2[1:]) * np.cos(np.arange(1, FDn) * kk)
    """
    % Find correlation bw/ Ecut (Ry) and mesh size h (Bohr)
    clear all; close all; clc;
    
    % input variables
    FDn = 6; % FD order / 2
    Ecut = 30; % in Ry
    
    
    
    
    
    
    epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.
    % Ecut in Ha
    %Ecut = Ecut * 0.5;
    
    % finite difference weights
    w2 = zeros(1,FDn+1); 
    for k=1:FDn
        w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
                        (k*k*factorial(FDn-k)*factorial(FDn+k));
        w2(1) = w2(1)-2*(1/(k*k));
    end
    
    kk = linspace(0,pi,1001);
    y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
    freq_err = abs(y_cos - kk.^2);
    kc = kk(max(find(freq_err < epsilon)))
    
    plot(kk,y_cos,'-',kk,kk.^2,'--','LineWidth',2);
    hold on
    plot(kc,-w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kc),'p','MarkerSize',15);
    hold off
    
    kc_by_pi = kc / pi;
    
    h = kc / sqrt(Ecut);
    
    fprintf('******************************************\n');
    fprintf('Ecut = %.2f Ry = %.2f Ha corresponds to:\n\t h = %.3f Bohr in real space\n', ...
            Ecut,0.5*Ecut,h);
    fprintf('******************************************\n');
    
    
    h_pw = pi / sqrt(Ecut)
    
    N_r = (h_pw / h)^3
    """


def h2Ecut():
    # to do, convert to python and test
    """
    % Find correlation bw/ Ecut (Ry) and mesh size h (Bohr)
    clear all; close all; clc;

    % input variables
    FDn = 6; % FD order / 2
    Ecut = 30; % in Ry






    epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.
    % Ecut in Ha
    %Ecut = Ecut * 0.5;

    % finite difference weights
    w2 = zeros(1,FDn+1); 
    for k=1:FDn
        w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
                        (k*k*factorial(FDn-k)*factorial(FDn+k));
        w2(1) = w2(1)-2*(1/(k*k));
    end

    kk = linspace(0,pi,1001);
    y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
    freq_err = abs(y_cos - kk.^2);
    kc = kk(max(find(freq_err < epsilon)))

    plot(kk,y_cos,'-',kk,kk.^2,'--','LineWidth',2);
    hold on
    plot(kc,-w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kc),'p','MarkerSize',15);
    hold off

    kc_by_pi = kc / pi;

    h = kc / sqrt(Ecut);

    fprintf('******************************************\n');
    fprintf('Ecut = %.2f Ry = %.2f Ha corresponds to:\n\t h = %.3f Bohr in real space\n', ...
            Ecut,0.5*Ecut,h);
    fprintf('******************************************\n');


    h_pw = pi / sqrt(Ecut)

    N_r = (h_pw / h)^3
    """


def cutoff2gridspacing(E):
    """Convert planewave energy cutoff to a real-space gridspacing.
       Taken from GPAW
    """
    return np.pi / np.sqrt(2 * E / Hartree) * Bohr
