from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator
import numpy as np
from ase import Atoms
from pymatgen.io.cif import CifWriter
from ase.io import read, write
from ase.build import surface, make_supercell
from ase.neighborlist import NeighborList
from ase.build import rotate
from pymatgen.analysis.local_env import CrystalNN
import os


def expand_slab(structure, expansion_matrix):
    """Expand the unit cell"""
    print("Before expansion:", len(structure))
    output = make_supercell(structure, expansion_matrix)
    output.write("expanded_structures.cif")
    return make_supercell(structure, expansion_matrix)
    #return structure * expansion_matrix


def reorient_surface(structure, miller_index, layers=1, vacuum=0):
    """Cut the surface according to the Miller index and add a vacuum layer"""
    lattice = Lattice(structure.get_cell())  # Get lattice matrix from ASE
    print("print lattice ********************")
    print(lattice)
    print("print lattice ********************")
    species = [atom.symbol for atom in structure]  # Get atomic species
    positions = structure.get_positions()  # Get atomic coordinates
    
    # Convert ASE Cartesian coordinates to pymatgen fractional coordinates
    frac_coords = lattice.get_fractional_coords(positions)  # Convert Cartesian to fractional coordinates
    #print(frac_coords) 
    #pymatgen_structure = Structure(lattice, species, positions)  # Create pymatgen structure object
    pymatgen_structure = Structure(lattice, species, frac_coords)  # Create pymatgen structure object
    angles = pymatgen_structure.lattice.angles  # Return alpha, beta, gamma
    print("Cell angles (degrees):", angles)
    latticepy = pymatgen_structure.lattice  # Get lattice matrix from pymatgen
    print(f"lattice information: {latticepy}")
    print('print convert structure *****')
    print(pymatgen_structure)

    # Output pymatgen structure as CIF file
    writer = CifWriter(pymatgen_structure)
    cif_filename = "converted_structure.cif"  # Adjust filename if needed
    writer.write_file(cif_filename)  # Save as CIF
    print(f"Converted structure saved as CIF: {cif_filename}")
    
    print(pymatgen_structure.lattice) 
    # Use pymatgen's SlabGenerator for cutting, ref https://pymatgen.org/pymatgen.core.html#pymatgen.core.surface.SlabGenerator
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, center_slab=True, lll_reduce=False, reorient_lattice=False)
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, lll_reduce=False, center_slab=False, reorient_lattice=False, primitive=False)
    slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, reorient_lattice=True, center_slab=False, primitive=False, max_normal_search=None)
    # print(slabgen)
    slab = slabgen.get_slab()
    #print(slabs)

    slab_cif_filename = "reorient_structure.cif"
    slab_writer = CifWriter(slab)
    slab_writer.write_file(slab_cif_filename)

    # Return the reoriented slab
    return slab

def cut_surface_with_constraints(structure, miller_index, layers=1, vacuum=20):
    """Cut the surface according to the Miller index and add a vacuum layer"""
    lattice = Lattice(structure.get_cell())  # Get lattice matrix from ASE
    print("print lattice ********************")
    print(lattice)
    print("print lattice ********************")
    species = [atom.symbol for atom in structure]  # Get atomic species
    positions = structure.get_positions()  # Get atomic coordinates
    
    # Convert ASE Cartesian coordinates to pymatgen fractional coordinates
    frac_coords = lattice.get_fractional_coords(positions)  # Convert Cartesian to fractional coordinates
    #print(frac_coords) 
    #pymatgen_structure = Structure(lattice, species, positions)  # Create pymatgen structure object
    pymatgen_structure = Structure(lattice, species, frac_coords)  # Create pymatgen structure object
    angles = pymatgen_structure.lattice.angles  # Return alpha, beta, gamma
    print("Cell angles (degrees):", angles)
    latticepy = pymatgen_structure.lattice  # Get lattice matrix from pymatgen
    print(f"lattice information: {latticepy}")
    print('print convert structure *****')
    print(pymatgen_structure)
    '''
    # Output pymatgen structure as CIF file
    writer = CifWriter(pymatgen_structure)
    cif_filename = "converted_structure.cif"  # Adjust filename if needed
    writer.write_file(cif_filename)  # Save as CIF
    print(f"Converted structure saved as CIF: {cif_filename}")
    '''    
    print(pymatgen_structure.lattice) 
    # Use pymatgen's SlabGenerator for cutting, ref https://pymatgen.org/pymatgen.core.html#pymatgen.core.surface.SlabGenerator
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, center_slab=True, lll_reduce=False, reorient_lattice=False)
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, lll_reduce=False, center_slab=False, reorient_lattice=False, primitive=False)
    slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, reorient_lattice=False, center_slab=False, primitive=False, max_normal_search=None)
    # print(slabgen)
    slabs = slabgen.get_slabs(tol=0.01,ftol=0.01)
    #print(slabs)
    for i, slab in enumerate(slabs):
        print(f"Processing slab {i+1}/{len(slabs)}")
    
        writer = CifWriter(slab)
        cif_filename = f"slabs_{i+1}.cif"  # Adjust filename if needed
        writer.write_file(cif_filename)  # Save as CIF

    # Return all generated slabs
    return slabs



def process_structure(cif_file, miller_index, expansion_matrix, layers=1, vacuum=20):
    """Main workflow: cut crystal surfaces, calculate bond counts, check termination validity, and save results"""
    structure = read(cif_file)
        
    # Print unit cell information
    print("Original cell parameters (a, b, c, alpha, beta, gamma):")
    print(structure.get_cell_lengths_and_angles())
        
    # Print number of atoms and species
    print(f"Number of atoms: {len(structure)}")
    print(f"Atomic species: {set(atom.symbol for atom in structure)}")
        
    # Print space group information (if available)
    spacegroup = structure.info.get('spacegroup', 'Not specified')
    print(f"Spacegroup: {spacegroup}")
    
    lattice = Lattice(structure.get_cell())  # Get lattice matrix from ASE
    print(f"lattice information: {lattice}")

    reorient_surface(structure, miller_index)
    
    structure = read('reorient_structure.cif')
    # Expand the supercell
    expanded_structure = expand_slab(structure, expansion_matrix)



    #structure = read('expanded_structures.cif')
    # Cut the surface and obtain all slabs
    slabs = cut_surface_with_constraints(
        expanded_structure, [0, 0, 1], layers=layers, vacuum=vacuum
        #structure, [0, 0, 1], layers=layers, vacuum=vacuum
    )
    '''
    # Iterate over all slabs
    for i, slab in enumerate(slabs):
        print(f"Processing slab {i+1}/{len(slabs)}")
    
        writer = CifWriter(slab)
        cif_filename = f"slabs_{miller_index}_{i+1}.cif"  # Adjust filename if needed
        writer.write_file(cif_filename)  # Save as CIF
    '''



# Example call
cif_file = "./Microcline-MS-1979.cif"  # Replace with actual CIF file path
miller_index = [1, 1, 0]  # Replace with actual Miller index
expansion_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]  # Replace with actual expansion matrix

process_structure(cif_file, miller_index, expansion_matrix)
