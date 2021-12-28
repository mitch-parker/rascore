# rascore: A tool for analyzing the conformations of RAS structures

## Summary

rascore is a Python package for classifying RAS (KRAS, NRAS, and HRAS) structures by their switch 1 (SW1; residues 25-40) and switch 2 (SW2: residues 56-76) loop conformations and annotating their molecular contents (e.g., mutation status, bound proteins and inhibitors). Details of our RAS conformational classification are provided on BioArxiv in the manuscript: "An expanded classification of active, inactive and druggable RAS conformations." We hope that researchers will use rascore to gain novel insights into RAS biology and drug discovery. 

## Installation

The following package versions are required:

- Bio==1.3.2
- rdkit==2009.Q1-1
- pandas==1.2.4
- numpy==1.20.1
- scikit_learn==1.0.1
- statsmodels==0.12.2
- matplotlib==3.3.4
- seaborn==0.11.1
- tqdm==4.59.0
- requests==2.25.1

Quickstart commands with an installation of Anaconda (https://www.anaconda.com/products/individual):

- conda install -c conda-forge -c schrodinger pymol-bundle
- conda install -c conda-forge rdkit
- conda install -c conda-forge cairosvg 
- conda install -c conda-forge tqdm
- conda install -c conda-forge statannot 
- conda install -c conda-forge matplotlib-venn

## Usage

Users may provide stuctures in mmCIF or PDB format, or rascore can be used to download available RAS structures from the Protein Data Bank (PDB). Importantly, any user inputted structure must be numbered according to UniProt scheme. Below, we provide commands for running the scripts in rascore:

# Authors

**Mitchell Parker**

- Email: mitch.isaac.parker@gmail.com
- GitHub: https://github.com/mitch-parker

**Roland Dunbrack**

- Email: roland.dunbrack@gmail.com
- GitHub: https://github.com/DunbrackLab

## License
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
