# -*- coding: utf-8 -*-
"""
  Copyright 2022 Mitchell Isaac Parker

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

import streamlit as st
from PIL import Image

from ..functions.table import mask_equal
from ..functions.col import pdb_code_col
from ..functions.path import pages_str, data_str, get_file_path
from ..functions.gui import load_st_table, write_st_end, create_st_button, show_st_structure, get_neighbor_path


def home_page():

    left_col, right_col = st.columns(2)

    df = load_st_table(__file__)

    show_st_structure(mask_equal(df, pdb_code_col, "6oim"),
            zoom=1.2,
            width=400,
            height=300,
            surface_trans=1,
            spin_on=True,
            st_col=left_col)

    right_col.markdown("# Rascore")
    right_col.markdown("### A tool for analyzing RAS protein structures")
    right_col.markdown("**Created by Mitchell Parker and Roland Dunbrack**")
    right_col.markdown("**Fox Chase Cancer Center**")

    database_link_dict = {
        "Cancer Research Paper": "https://aacrjournals.org/cancerres/article/doi/10.1158/0008-5472.CAN-22-0804/696349/Delineating-The-RAS-Conformational-LandscapeThe",
        "GitHub Page": "https://github.com/mitch-parker/rascore",
        "RCSB Protein Data Bank": "https://www.rcsb.org",
    }

    st.sidebar.markdown("## Database-Related Links")
    for link_text, link_url in database_link_dict.items():
        create_st_button(link_text, link_url, st_col=st.sidebar)

    community_link_dict = {
        "NCI RAS Initiative": "https://www.cancer.gov/research/key-initiatives/ras",
        "KRAS Kickers": "https://www.kraskickers.org",
        "RASopathies Network": "https://rasopathiesnet.org",
    }

    st.sidebar.markdown("## Community-Related Links")
    for link_text, link_url in community_link_dict.items():
        create_st_button(link_text, link_url, st_col=st.sidebar)

    software_link_dict = {
        "BioPython": "https://biopython.org",
        "RDKit": "https://www.rdkit.org",
        "PDBrenum": "http://dunbrack.fccc.edu/PDBrenum/",
        "Fpocket": "https://bioserv.rpbs.univ-paris-diderot.fr/services/fpocket/",
        "PyMOL": "https://pymol.org/2/",
        "3Dmol": "https://3dmol.csb.pitt.edu",
        "Pandas": "https://pandas.pydata.org",
        "NumPy": "https://numpy.org",
        "SciPy": "https://scipy.org",
        "Sklearn": "https://scikit-learn.org/stable/",
        "Matplotlib": "https://matplotlib.org",
        "Seaborn": "https://seaborn.pydata.org",
        "Streamlit": "https://streamlit.io",
    }

    st.sidebar.markdown("## Software-Related Links")
    link_1_col, link_2_col, link_3_col = st.sidebar.columns(3)

    i = 0
    link_col_dict = {0: link_1_col, 1: link_2_col, 2: link_3_col}
    for link_text, link_url in software_link_dict.items():

        st_col = link_col_dict[i]
        i += 1
        if i == len(link_col_dict.keys()):
            i = 0

        create_st_button(link_text, link_url, st_col=st_col)

    st.markdown("---")

    st.markdown(
        """
        ### Summary
        *Rascore* is a tool for analyzing structures of the RAS protein family
        (KRAS, NRAS, and HRAS). The *Rascore* 
        database presents a continually updated analysis of all available
        RAS structures in the Protein Data Bank (PDB) with their catalytic switch 1 (SW1) 
        and switch 2 (SW2) loops conformationally classified and their molecular 
        contents annotated (e.g., mutation status, nucleotide state, 
        bound protein, inhibitor site). 

        Details of our work are 
        provided in the [*Cancer Research*](https://aacrjournals.org/cancerres/article/doi/10.1158/0008-5472.CAN-22-0804/696349/Delineating-The-RAS-Conformational-LandscapeThe)
        paper, **Delineating The RAS Conformational Landscape**.
        We hope that researchers will use 
        *Rascore* to gain novel insights into RAS biology and 
        drug discovery. 
        """
    )

    left_col, right_col = st.columns(2)

    img = Image.open(
        get_file_path(
            "rascore_abstract.png",
            dir_path=get_neighbor_path(__file__, pages_str, data_str),
        )
    )

    right_col.image(img, output_format="PNG")


    left_col.markdown(
        """
        ### Usage

        To the left, is a dropdown main menu for navigating to 
        each page in the *Rascore* database:

        - **Home Page:** We are here!
        - **Database Overview:** Overview of the *Rascore* database, molecular annotations, and RAS conformational classification.
        - **Search PDB:** Search for individual PDB entries containing RAS structures.
        - **Explore Conformations:** Explore RAS SW1 and SW2 conformations found in the PDB by nucleotide state.
        - **Analyze Mutations:** Analyze the structural impact of RAS mutations by comparing WT and mutated structures.
        - **Compare Inhibitors:** Compare inhibitor-bound RAS structures by compound binding site and chemical substructure.
        - **Query Database:** Query the *Rascore* database by conformations and molecular annotations.
        - **Classify Structures:** Conformationally classify and annotate the molecular contents of uploaded RAS structures.
        """
    )
    st.markdown("---")

    left_info_col, right_info_col = st.columns(2)

    left_info_col.markdown(
        f"""
        ### Authors
        Please feel free to contact us with any issues, comments, or questions.

        ##### Mitchell Parker [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40Mitch_P)](https://twitter.com/Mitch_P)

        - Email:  <mip34@drexel.edu> or <mitchell.parker@fccc.edu>
        - GitHub: https://github.com/mitch-parker

        ##### Roland Dunbrack [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40RolandDunbrack)](https://twitter.com/RolandDunbrack)

        - Email: <roland.dunbrack@fccc.edu>
        - GitHub: https://github.com/DunbrackLab
        """,
        unsafe_allow_html=True,
    )

    right_info_col.markdown(
        """
        ### Funding

        - NIH NIGMS F30 GM142263 (to M.P.)
        - NIH NIGMS R35 GM122517 (to R.D.)
         """
    )

    right_info_col.markdown(
        """
        ### License
        Apache License 2.0
        """
    )

    write_st_end()
