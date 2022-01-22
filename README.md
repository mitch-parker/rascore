![alt text](https://github.com/mitch-parker/rascore/blob/main/rascore/data/rascore_logo.png?)

# Summary

*rascore* is a tool for analyzing RAS structures (KRAS, NRAS, and HRAS) by the conformations of their catalytic switch 1 (SW1) and switch 2 (SW2) loops. In addition, *rascore* can be used to build and query an updatable database of all available RAS structures from the Protein Data Bank with their SW1 and SW2 loops conformationally classified and their molecular contents annotated (e.g., mutation status, nucleotide state, bound proteins, small molecule inhibitors, etc.). 

Details of our RAS conformational classification and approach are provided on *BioArxiv* in the manuscript: **An expanded classification of active, inactive and druggable RAS conformations.** We hope that researchers will use *rascore* to gain novel insights into RAS biology and drug discovery. 

# Installation

**Quickstart environment setup with an installation of Anaconda (https://www.anaconda.com/products/individual):**

pip install rascore

conda create -n rascore_env

conda activate rascore_env

conda install -n rascore_env -c conda-forge -c schrodinger pymol-bundle

conda install -c conda-forge rdkit

conda install -c conda-forge cairosvg

conda install -c conda-forge tqdm

conda install -c conda-forge statannot

conda install -c conda-forge matplotlib-venn

conda install -c conda-forge fpocket

conda install -c conda-forge pyfiglet

pip install streamlit

pip install futures

pip install py3Dmol

pip install stmol

# Usage

### 1) Conformationally classify RAS structure(s):

rascore -classify *[path to coordinate file(s)]* -out *[output directory path]*

*Note.* Coordinate files must be provided in mmCIF or PDB formats and numbered according to UniProt scheme. The following inputs are possible: 
- Space-separated list
- Line-separated list file
- Tab-separated table file with columns core_path (coordinate path), modelid (optional, model number), chainid (chain identifier), nuc_class (optional, nucleotide class)

### 2) Build rascore database from Protein Data Bank:

rascore -build *[optional, path to pdbaa file]* -out *[output directory path]*

*Note.* Take ~1 hour to build from scratch and requires ~3 GB of storage.

### 3) Run rascore GUI application:

rascore -app *[optional, path to rascore database directory]* -out *[output directory path]*

*Note.* Can run limited version if rascore database directory not specified.

# Authors

Please feel free to contact us with any issues, comments, or questions regarding rascore.

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

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.