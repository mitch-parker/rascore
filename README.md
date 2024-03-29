![alt text](https://github.com/mitch-parker/rascore/blob/main/src/rascore/util/data/rascore_logo.png)

<a href="https://github.com/mitch-parker/rascore" title="Go to GitHub repo"><img src="https://img.shields.io/static/v1?label=mitch-parker&message=rascore&color=e78ac3&logo=github" alt="mitch-parker - rascore"></a>
<a href="https://github.com/mitch-parker/rascore"><img src="https://img.shields.io/github/stars/mitch-parker/rascore?style=social" alt="stars - rascore"></a>
<a href="https://github.com/mitch-parker/rascore"><img src="https://img.shields.io/github/forks/mitch-parker/rascore?style=social" alt="forks - rascore"></a>

</div>

<a href="https://github.com/mitch-parker/rascore/releases/"><img src="https://img.shields.io/github/release/mitch-parker/rascore?include_prereleases=&sort=semver&color=e78ac3" alt="GitHub release"></a>
<a href="#license"><img src="https://img.shields.io/badge/License-Apache_2.0-e78ac3" alt="License"></a>
<a href="https://github.com/mitch-parker/rascore/issues"><img src="https://img.shields.io/github/issues/mitch-parker/rascore" alt="issues - rascore"></a>

## Summary

*Rascore* is a tool for analyzing structures of the RAS protein family (KRAS, NRAS, and HRAS). In addition, *Rascore* can be used to build an updatable database of all available RAS structures in the Protein Data Bank (PDB) with their catalytic switch 1 (SW1) and switch 2 (SW2) loops conformationally classified and their 
molecular contents annotated (e.g., mutation status, nucleotide state, bound protein, inhibitor site). 

Details of our work are provided in the [*Cancer Research*](https://aacrjournals.org/cancerres/article/doi/10.1158/0008-5472.CAN-22-0804/696349/Delineating-The-RAS-Conformational-LandscapeThe) paper, **Delineating The RAS Conformational Landscape**. We hope that researchers will use *Rascore* to gain novel insights into RAS biology and drug discovery. 

## Web Database

A continually updated version of the *Rascore* database is hosted at http://dunbrack.fccc.edu/rascore/ or https://rascore.streamlitapp.com.

## Local Installation

**Quickstart commands for environment setup with [Anaconda](https://www.anaconda.com/products/individual):**

```conda create -n rascore_env python=3.8```

```conda activate rascore_env```

```pip install rascore```

```conda install -c conda-forge -c schrodinger pymol-bundle=2.5.2```

```conda install -c conda-forge fpocket=4.0.0```

(Additional for macOS: pip install PyQt5)

## Usage

*Note.* Must activate Anaconda enironment with command "conda activate rascore_env" before running.

#### 1. Conformationally classify RAS structures:

```rascore -classify [path to coordinate files] -out [output directory path]```

- Coordinate files must be provided in mmCIF or PDB formats and numbered according to UniProt scheme. The following inputs are possible: 
    - Space-separated list
    - Line-separated list file
    - Tab-separated table file with columns core_path (coordinate path), modelid (optional, model number), chainid (chain identifier), nucleotide_class (optional, nucleotide state)

#### 2. Build the Rascore database from the PDB:

```rascore -build [optional, path to pdbaa file] -out [output directory path]```

- Takes 1 hour to build from scratch and requires 3 GB of storage.

#### 3. Run the Rascore database application locally:

```rascore -gui [optional, path to the Rascore database directory] -out [output directory path]```

## Authors

Please feel free to contact us with any issues, comments, or questions.

#### Mitchell Parker [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40Mitch_P)](https://twitter.com/Mitch_P)

- Email: <mip34@drexel.edu> or <mitchell.parker@fccc.edu>
- GitHub: https://github.com/mitch-parker

#### Roland Dunbrack [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40RolandDunbrack)](https://twitter.com/RolandDunbrack)

- Email: <roland.dunbrack@fccc.edu>
- GitHub: https://github.com/DunbrackLab

## Funding

- NIH NIGMS F30 GM142263 (to M.P.)
- NIH NIGMS R35 GM122517 (to R.D.)

## License
Apache License 2.0


Copyright (c) 2022 Mitchell Isaac Parker