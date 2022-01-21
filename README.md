# rascore: A tool for analyzing the conformations of RAS structures

## Summary

rascore is a package for conformationally classifying RAS structures (KRAS, NRAS, and HRAS) by their switch 1 and switch 2 loop conformations. Details of our RAS conformational classification are provided on BioArxiv in the manuscript: "An expanded classification of active, inactive and druggable RAS conformations." We hope that researchers will use rascore to gain novel insights into RAS biology and drug discovery. 

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
- conda install -c conda-forge fpocket
- conda install -c conda-forge pyfiglet 

## Usage

Users may provide RAS stuctures in mmCIF or PDB format for conformational classification, or rascore can be used to download all available KRAS, NRAS, and HRAS structures from the Protein Data Bank (PDB) to conformationally classify. Importantly, any user inputted structure(s) must be numbered according to UniProt scheme. Below, we provide commands for running the scripts in rascore:

# Authors

**Mitchell Parker**

- Email: mitch.isaac.parker@gmail.com
- GitHub: https://github.com/mitch-parker

**Roland Dunbrack**

- Email: roland.dunbrack@gmail.com
- GitHub: https://github.com/DunbrackLab

## License
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