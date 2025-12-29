# Data from "Molecular mechanism of heterogeneous ice nucleation in the atmosphere"
## Authors: Wanqi Zhou and Pablo M. Piaggi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17610333.svg)](https://doi.org/10.5281/zenodo.17610332)

—————————————————————————————————————————

This repository contains the scripts for searching slab surfaces, the input file for density functional theory (DFT) calculations used to generate the training dataset, input files for training machine-learning potentials (MLPs) and running molecular dynamics (MD) simulations. It also includes the generated datasets and the final trained MLP models.

In this README, we introduce each file according to our four-step workflow: (1) slab surface generation, (2) DFT calculations, (3) MLP training, and (4) MD simulations.

## (1) Slab surface generation
All possible terminations of the (001), (010), (100), (110), ($\bar{1}$10), and ($\bar{2}$01) planes were generated using the \textsc{SlabGenerator} module of the \textsc{pymatgen} package\cite{tran2016surface}.
A cutoff of 1~\text{\AA} along the surface normal direction was applied to determine whether atoms lie within the same atomic plane, which was used to filter the generated terminations.


This repository contains the scripts, input files, and workflow used to generate the slab surfaces, prepare the training dataset, train the machine-learning potentials (MLPs), and run the molecular dynamics (MD) and density functional theory (DFT) calculations reported in our work.


📦 Dataset & MLP Download

Due to file size, the training dataset and the final trained MLP model are stored on Google Drive.

📁 Training Dataset

Contains raw atomic configurations and DFT energies/forces for MLP training.

👉 Download: https://drive.google.com/drive/folders/1qKnz3tHYAP0c35sSq0amDg7CQuScfyqp?usp=drive_link

🤖 Final Trained MLP

The full trained MLP model, ready for MD simulations.

👉 Download:https://drive.google.com/drive/folders/1iJiQLxTOqFKbddP4vH3R3-Vdb4NTtjK8?usp=drive_link

