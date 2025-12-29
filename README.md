# Data from "Molecular mechanism of heterogeneous ice nucleation in the atmosphere"
## Authors: Wanqi Zhou and Pablo M. Piaggi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17610333.svg)](https://doi.org/10.5281/zenodo.17610332)

—————————————————————————————————————————

This repository contains the scripts for searching slab surfaces, the input file for density functional theory (DFT) calculations used to generate the training dataset, input files for training machine-learning potentials (MLPs) and running molecular dynamics (MD) simulations. It also includes the generated datasets and the final trained MLP models.

In this README, we introduce each file according to our four-step workflow: (1) slab surface generation, (2) DFT calculations, (3) MLP training, and (4) MD simulations.

## (1) Slab surface generation
All possible terminations of the (001), (010), (100), (110), (−110), and (−201) crystallographic planes were generated using the SlabGenerator module from the pymatgen package.
The script Rotate_expand_cut_pymatgen.py is located at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Generating_slab_surfaces
This script performs the following steps:

1) Read crystal structure: The script first reads the CIF file Microcline-MS-1979.cif, which contains the feldspar unit cell, and prints basic information about the crystal structure. (output file, converted_structure.cif: structure converted from ASE to pymatgen format)

2) Rotate the unit cell: The crystal cell is rotated so that the specified Miller index is aligned with the z-axis. (output file, reorient_structure.cif: structure after Miller-index reorientation)

3) Expand the supercell (optional): The unit cell can be expanded along selected directions to construct a supercell, depending on the simulation requirements. (output file, expanded_structures.cif: expanded supercell)

4) Cut the slab: Finally, the structure is cut along the z-axis to generate the desired slab geometry.
A cutoff of 1 Å along the surface-normal direction was applied to determine whether atoms belong to the same atomic plane.
This criterion was used to filter and identify unique slab terminations.
The script outputs all possible slab surface terminations that satisfy these conditions. (output file, slabs_*.cif: all generated slab surface terminations)

After running the Python script above, hydrogen atoms were added near undercoordinated oxygen atoms at the surfaces.
Additionally, surface terminations that disrupted the integrity of SiO₄ and AlO₄ polyhedra or violated charge neutrality were discarded.





This repository contains the scripts, input files, and workflow used to generate the slab surfaces, prepare the training dataset, train the machine-learning potentials (MLPs), and run the molecular dynamics (MD) and density functional theory (DFT) calculations reported in our work.


📦 Dataset & MLP Download

Due to file size, the training dataset and the final trained MLP model are stored on Google Drive.

📁 Training Dataset

Contains raw atomic configurations and DFT energies/forces for MLP training.

👉 Download: https://drive.google.com/drive/folders/1qKnz3tHYAP0c35sSq0amDg7CQuScfyqp?usp=drive_link

🤖 Final Trained MLP

The full trained MLP model, ready for MD simulations.

👉 Download:https://drive.google.com/drive/folders/1iJiQLxTOqFKbddP4vH3R3-Vdb4NTtjK8?usp=drive_link

