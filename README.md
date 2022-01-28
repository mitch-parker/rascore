![alt text](https://github.com/mitch-parker/rascore/blob/main/src/rascore/util/data/rascore_logo.png)

<a href="https://github.com/mitch-parker/rascore" title="Go to GitHub repo"><img src="https://img.shields.io/static/v1?label=mitch-parker&message=rascore&color=%23e78ac3&logo=github" alt="mitch-parker - rascore"></a>
<a href="https://github.com/mitch-parker/rascore"><img src="https://img.shields.io/github/stars/mitch-parker/rascore?style=social" alt="stars - rascore"></a>
<a href="https://github.com/mitch-parker/rascore"><img src="https://img.shields.io/github/forks/mitch-parker/rascore?style=social" alt="forks - rascore"></a>
</div>
<a href="https://github.com/mitch-parker/rascore/releases/"><img src="https://img.shields.io/github/release/mitch-parker/rascore?include_prereleases=&sort=semver&color=e78ac3" alt="GitHub release"></a>
<a href="#license"><img src="https://img.shields.io/badge/License-MIT-e78ac3" alt="License"></a>


# Summary

*rascore* is a tool for analyzing RAS structures (KRAS, NRAS, and HRAS) by the conformations of their catalytic switch 1 (SW1) and switch 2 (SW2) loops. In addition, *rascore* can be used to search and query an updatable database of all available RAS structures in the Protein Data Bank with their SW1 and SW2 loops conformationally classified and their molecular contents annotated (e.g., mutation status, nucleotide state, bound proteins, small molecule inhibitors, etc.). 

Details of our RAS conformational classification are provided on *BioArxiv* in the manuscript: **An expanded classification of active, inactive and druggable RAS conformations.** We hope that researchers will use *rascore* to gain novel insights into RAS biology and drug discovery. 

# Graphical User Interface (GUI) Application

A continually updated version of the *rascore* database is hosted on Streamlit at https://share.streamlit.io/mitch-parker/rascore/main/src/rascore/rascore_gui.py. A version of the GUI applicaiton with extra user functionality can be run with a local installation.

# Local Installation

**Quickstart commands for environment setup with [Anaconda](https://www.anaconda.com/products/individual):**

conda create -n rascore_env python=3.8

conda activate rascore_env

pip install rascore 

conda install -c conda-forge -c schrodinger pymol-bundle=2.5.2

conda install -c conda-forge fpocket=4.0.0

(Additional for macOS: pip install PyQt5)

# Usage

*Note.* Must activate Anaconda enironment with command "conda activate rascore_env" before running.

### 1) Conformationally classify RAS structure(s):

rascore -classify *[path to coordinate file(s)]* -out *[output directory path]*

- Coordinate files must be provided in mmCIF or PDB formats and numbered according to UniProt scheme. The following inputs are possible: 
    - Space-separated list
    - Line-separated list file
    - Tab-separated table file with columns core_path (coordinate path), modelid (optional, model number), chainid (chain identifier), nuc_class (optional, nucleotide class)

### 2) Build the rascore database from the Protein Data Bank:

rascore -build *[optional, path to pdbaa file]* -out *[output directory path]*

- Take ~1 hour to build from scratch and requires ~3 GB of storage.

### 3) Run the rascore GUI application:

rascore -gui *[optional, path to the rascore database directory]* -out *[output directory path]*

- Will include extra user functionality if the rascore database directory path is specified.

# Authors

Please feel free to contact us with any issues, comments, or questions.

### Mitchell Parker [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40Mitch_P)](https://twitter.com/Mitch_P)

- Email: <mitch.isaac.parker@gmail.com>
- GitHub: https://github.com/mitch-parker

### Roland Dunbrack [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40RolandDunbrack)](https://twitter.com/RolandDunbrack)

- Email: <roland.dunbrack@gmail.com>
- GitHub: https://github.com/DunbrackLab

# Funding

- NIH F30 GM142263 (to M.P.)
- NIH R35 GM122517 (to R.D.)

# License
MIT License

Copyright (c) 2022 Mitchell Isaac Parker