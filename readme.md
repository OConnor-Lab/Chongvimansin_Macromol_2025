Melt Scission Analysis

This repository contains the scripts and data necessary for analyzing polymer melt scission simulations and generating plots for the associated manuscript.

Repository Structure

meltscission.py: A Python script used for processing simulation data and extracting relevant scission information.

Analysis_Plots.ipynb: A Jupyter Notebook that provides a guided analysis workflow, including data processing and visualization.

data/: A directory containing the processed data files required for generating plots. Note: Raw simulation data files are not included due to space constraints.

LAMMPS Version and Required Packages

This analysis is based on LAMMPS version lammps-29Sep2021. The following LAMMPS packages are required:

MOLECULE

UEF

UEFEX

NETCDF

Additionally, modified bond and angle .cpp and .h files that allow bond scission are uploaded in this repository.

Requirements

The scripts require the following Python packages:

numpy

matplotlib

scipy

jupyter

Install dependencies using:

pip install numpy pandas matplotlib scipy jupyter

Usage Guide

Prepare Data: Ensure the processed data files are placed inside the data/ folder.

Run the Notebook: Open and execute Analysis_Plots.ipynb to generate the figures for the manuscript.

Standalone Script Execution: Run meltscission.py to process new simulation data before using the notebook.

python meltscission.py

Notes

The Jupyter notebook provides step-by-step guidance on the analysis process.

Modify paths in the notebook if necessary to match your directory structure.

If additional simulation data needs to be processed, update meltscission.py accordingly.

Citation

If you use this code or data in your research, please cite the following manuscript:



For any questions or issues, please contact the project maintainer.