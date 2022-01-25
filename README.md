![alt text](https://github.com/mitch-parker/rascore/blob/main/src/rascore/scripts/data/rascore_logo.png?)

# Summary

*rascore* is a tool for analyzing RAS structures (KRAS, NRAS, and HRAS) by the conformations of their catalytic switch 1 (SW1) and switch 2 (SW2) loops. In addition, *rascore* can be used to build and query an updatable database of all available RAS structures from the Protein Data Bank with their SW1 and SW2 loops conformationally classified and their molecular contents annotated (e.g., mutation status, nucleotide state, bound proteins, small molecule inhibitors, etc.). 

Details of our RAS conformational classification are provided on *BioArxiv* in the manuscript: **An expanded classification of active, inactive and druggable RAS conformations.** We hope that researchers will use *rascore* to gain novel insights into RAS biology and drug discovery. 

# Installation

**Quickstart environment setup with an installation of Anaconda (https://www.anaconda.com/products/individual):**

conda create -n rascore_env python=3.8

conda activate rascore_env

pip install rascore==0.1.25

conda install -c conda-forge -c schrodinger pymol-bundle=2.5.2

conda install -c conda-forge fpocket=4.0.0

# Usage

### 1) Conformationally classify RAS structure(s):

rascore -classify *[path to coordinate file(s)]* -out *[output directory path]*

*Note.* Coordinate files must be provided in mmCIF or PDB formats and numbered according to UniProt scheme. The following inputs are possible: 
- Space-separated list
- Line-separated list file
- Tab-separated table file with columns core_path (coordinate path), modelid (optional, model number), chainid (chain identifier), nuc_class (optional, nucleotide class)

### 2) Build the rascore database from the Protein Data Bank:

rascore -build *[optional, path to pdbaa file]* -out *[output directory path]*

*Note.* Take ~1 hour to build from scratch and requires ~3 GB of storage.

### 3) Run the rascore GUI application (IN PROGRESS):

rascore -gui *[optional, path to the rascore database directory]* -out *[output directory path]*

*Note.* Can run limited version if the rascore database directory path is not specified.

# Authors

Please feel free to contact us with any issues, comments, or questions.

### Mitchell Parker

- Email: <mitch.isaac.parker@gmail.com>
- GitHub: https://github.com/mitch-parker

### Roland Dunbrack

- Email: <roland.dunbrack@gmail.com>
- GitHub: https://github.com/DunbrackLab

# Funding

- NIH F30 GM142263 (to M.P.)
- NIH R35 GM122517 (to R.D.)

# License
MIT License

Copyright (c) 2022 Mitchell Isaac Parker