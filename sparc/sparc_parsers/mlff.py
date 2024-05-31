"""
Created on Fri May 24th 2024

Lucas Timmerman (Georgia Tech)

This file is provides the core functionality for interacting with the SPARC MLFF routines.
"""
from warnings import warn

import numpy as np
from ase.units import AUT, Angstrom, Bohr, GPa, Hartree, fs

def write_structures_file(structures:list, filename:str = "MLFF_data_reference_structures.txt"):
    """
    Write the structures to a file in the format that the MLFF code expects.
    """
    assert len(structures) > 0, "No structures provided to write to file."
    
    with open(filename, "w") as f:
        f.write("SPARC MLFF on-the-fly training data\n")
        f.write("Structures:\n")
    
        for i, atoms in enumerate(structures):
            f.write(f"structure_no: {i+1}\n")
            # Write cell information
            cell = atoms.get_cell() / Bohr  # Convert from Angstrom to Bohr
            mags = atoms.cell.lengths() / Bohr  # Convert from Angstrom to Bohr
            f.write("CELL:\n")
            f.write(f"{mags[0]} {mags[1]} {mags[2]}\n")
            
            # Write lattice unit vectors
            f.write("LatUVec:\n")
            for idx,vec in enumerate(cell):
                for _element in vec:
                    f.write(f"{_element/mags[idx]} ")
                f.write("\n")
            
            # Write number of atoms
            natoms = len(atoms)
            f.write(f"natom:\n{natoms}\n")
            
            # Write number of each type of atom
            symbols = atoms.get_chemical_symbols()
            unique_symbols = set(symbols)
            f.write("natom_elem:\n")
            for symbol in unique_symbols:
                f.write(f"{symbols.count(symbol)}\n")
            
            # Write atom positions
            positions = atoms.get_positions() / Bohr  # Convert from Angstrom to Bohr
            f.write("Atom Positions:\n")
            for pos in positions:
                f.write(f"{pos[0]} {pos[1]} {pos[2]}\n")
            
            # Write total energy
            energy = atoms.get_potential_energy() / Hartree  # Convert from eV to Hartree
            f.write(f"Etot(Ha)\n{energy}\n")
            
            # Write forces
            forces = atoms.get_forces() / (Hartree / Bohr)  # Convert from eV/Angstrom to Ha/Bohr
            f.write("F(Ha/Bohr)\n")
            for force in forces:
                f.write(f"{force[0]} {force[1]} {force[2]}\n")
            
            # Write stress if available
            if atoms.has('stress'):
                stress = atoms.get_stress(voigt=False) / GPa  # Convert from eV/Angstrom^3 to GPa
                f.write("Stress\n")
                for row in stress:
                    f.write(f"{row[0]} {row[1]} {row[2]}\n")
            
