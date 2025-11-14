# ice-nucleation-on-feldspar
Input files for manuscript
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
This repository contains the scripts, input files, and workflow used to generate the slab surfaces, prepare the training dataset, train the machine-learning potentials (MLPs), and run the molecular dynamics (MD) and density functional theory (DFT) calculations reported in our work.

.
├── scripts/                # Scripts for slab generation, preprocessing, analysis
├── training_input/         # Input files for training the MLP
├── validation_input/       # Validation data and scripts
├── md_simulations/         # Input files for MD simulations using the trained MLP
├── dft_calculations/       # Example DFT input files
└── README.md

📦 Dataset & MLP Download

Due to file size, the training dataset and the final trained MLP model are stored on Google Drive.

📁 Training Dataset

Contains raw atomic configurations, DFT energies/forces, and metadata used for MLP training.

👉 Download:
Dataset_for_training_MLP (Google Drive)
