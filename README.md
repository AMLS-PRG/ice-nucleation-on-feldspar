# Data from "Molecular mechanism of heterogeneous ice nucleation in the atmosphere"
## Authors: Wanqi Zhou and Pablo M. Piaggi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17610333.svg)](https://doi.org/10.5281/zenodo.17610332)

—————————————————————————————————————————

This repository contains the scripts for searching slab surfaces, the input file for density functional theory (DFT) calculations used to generate the training dataset, input files for training machine-learning potentials (MLPs) and running molecular dynamics (MD) simulations. It also includes the generated datasets and the final trained MLP models.

In this README, we introduce each file according to our four-step workflow: (1) slab surface generation, (2) DFT calculations, (3) MLP training, and (4) MD simulations.

## (1) 🏗️ Slab surface generation
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


## (2) 📦 DFT calculations (Dataset)
Plane-wave DFT calculations were performed using Quantum ESPRESSO v6.4.1.
An example input file, pw-water-0.in, is available at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Input_files_for_DFT

Key computational settings:

1) Exchange-correlation functional: SCAN (Strongly Constrained and Appropriately Normed) evaluated via LIBXC 4.3.4.

2) Pseudopotentials: Norm-conserving, scalar-relativistic for K, Al, Si, O, and H, parametrized using PBE, with 9, 11, 4, 6, and 1 valence electrons, respectively.

3) Kinetic energy cutoffs: 110 Ry for wave functions and 440 Ry for charge density.

4) k-point sampling: Only the Γ-point was used.

5) All other parameters were set to their default values in Quantum ESPRESSO.

Due to file size, the training dataset containing raw atomic configurations and DFT energies/forces for MLP training is stored on Google Drive. 👉 Download: https://drive.google.com/drive/folders/1qKnz3tHYAP0c35sSq0amDg7CQuScfyqp?usp=drive_link 📁

## (3) 🤖 MLP training (Final trained MLP)
The smooth edition of the Deep Potential methodology developed by Zhang et al., as implemented in DeePMD-kit v2.10.0, was used to train the machine-learning interatomic potentials (MLPs).
A typical input file, input.json, is available at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Input_files_for_training_MLP

The settings of the DeePMD-kit are as follows:

1) The sizes of the embedding and fitting networks were (50, 100, 200) and (120, 120, 120), respectively.

2) A smooth and hard cutoff radius of 6 Å and 3 Å were used.

3) The hyperparameters start_pref_e, start_pref_f, limit_pref_e, and limit_pref_f, which control the relative weights of energy and force terms in the total loss function, were set to 0.02, 1000, 1.0, and 10.0, respectively.

4) The initial learning rate was 0.002, with a decay step of 20,000, and the total number of training steps was 2 × 10⁶.

An active learning strategy was employed during the training process:
MD simulations driven by our previous MLP were performed to generate a series of configurations of the water-feldspar interfaces, covering all 13 K-feldspar terminations.
The energies and forces for these configurations were then calculated using SCAN DFT to further expand the dataset.
The resulting dataset, which included new configurations and their corresponding energies and atomic forces, was used to train a set of four MLPs.
Based on the newly trained MLPs, additional configurations were explored, and this cycle was repeated iteratively until a high-accuracy MLP with SCAN-level precision was obtained.

Due to file size, the final trained MLP is stored on Google Drive. 👉 Download: https://drive.google.com/drive/folders/1iJiQLxTOqFKbddP4vH3R3-Vdb4NTtjK8?usp=drive_link 📁

## (4) 🏃‍♂️ MD simulations
Molecular dynamics (MD) and enhanced sampling simulations were performed using LAMMPS interfaced with the PLUMED plugin.
All simulations employed the custom-trained MLPs described above to model atomic interactions.

### Standard MD simulations
An example input file for standard MD simulations is available at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Input_files_for_MD/Standard_MD

Simulation Details

1) Atomic masses: Standard atomic weights were used for K, Al, Si, and O. Hydrogen was assigned a mass of 2 g/mol to improve stability of the integration.

2) Boundary conditions: Periodic in all directions.

3) Time step: 0.5 fs.

4) Energy minimization: Conjugate gradient (CG) method.

5) Equilibration: NVT ensemble with a stochastic velocity-rescaling thermostat (0.1 ps relaxation time).

6) Simulation temperature: 290 K (supercooling ΔT = 18 K) for studies of water structure on multiple surfaces.

7) Simulation duration: More than 3 ns.

8) Constraints: K, Al, and Si atoms were constrained to their initial positions using a harmonic spring (force constant = 20 eV/Å²).

9) Trajectory output: Saved every 1000 steps.

### Steered MD simulations
As standard MD simulations cannot capture the nucleation process within affordable simulation times, an enhanced sampling method was employed: steered MD simulations guided by the Q₆ Steinhardt order parameter, as implemented in PLUMED (file: plumed.dat).
An example input file for steered MD simulations is available at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Input_files_for_MD/SteeredMD

Simulation Details

1) Order parameter: Calculated using the oxygen atoms of water molecules to quantify local structural ordering.

2) Restraint: A moving restraint with a force constant of 2 × 10⁶ kJ/mol applied to the collective variable.

3) Simulation temperature: 300 K (supercooling ΔT = 8 K).

4) Restraint protocol: Target value gradually increased from 0.05 to 0.2 over a 5 ns simulation.

5) Output:
Order parameter recorded every 1000 steps to monitor nucleation progress.
Trajectory saved every 10,000 steps.

### Umbrella sampling method
The free energy profile of ice nucleation was calculated using the umbrella sampling method. 
An example input file for umbrella sampling method is available at:
https://github.com/AMLS-PRG/ice-nucleation-on-feldspar/tree/main/Input_files_for_MD/US

Simulation Details

1) Collective variable: Q₆ Steinhardt parameter.

2) Windows: 41 windows in the range 0.05 ≤ Q₆ ≤ 0.11 with a width of 0.0015.

3) Restraint: Harmonic potential with a spring constant of 1 × 10⁶ kJ/mol applied to Q₆ in each window.

4) Simulation duration: 3 ns per window.

5) Free energy reconstruction: The time series of Q₆ from all windows were used to construct the free energy profile via the weighted histogram analysis method (WHAM).







