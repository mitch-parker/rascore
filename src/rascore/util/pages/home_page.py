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

from ..constants.pharm import pharm_color_dict, sp2_name
from ..constants.conf import loop_resid_dict, sw1_name, sw2_name, sw1_color, sw2_color
from ..functions.gui import write_st_end, create_st_button, show_st_structure


def home_page():

    style_lst = list()
    surface_lst = list()

    opacity = 1

    surface_lst = [
        [
            {"opacity": opacity, "color": "white"},
            {"chain": "A", "hetflag": False},
        ]
    ]

    style_lst.append(
        [
            {
                "resn": "GDP",
            },
            {
                "stick": {
                    "colorscheme": "whiteCarbon",
                    "radius": 0.2,
                }
            },
        ]
    )

    style_lst.append(
        [
            {
                "resn": "MG",
            },
            {
                "sphere": {
                    "color": "chartreuse",
                    "radius": 0.8,
                }
            },
        ]
    )

    style_lst.append(
        [
            {"chain": "A", "resi": 12, "atom": "CA"},
            {"sphere": {"color": "red", "radius": 0.8}},
        ]
    )

    style_lst.append(
        [
            {
                "chain": "A",
                "resn": "MOV",
                "elem": "C",
            },
            {"stick": {"color": pharm_color_dict[sp2_name], "radius": 0.2}},
        ]
    )
    style_lst.append(
        [
            {
                "chain": "A",
                "resn": "MOV",
                "elem": ["N", "O", "H"],
            },
            {"stick": {"colorscheme": "whiteCarbon", "radius": 0.2}},
        ]
    )

    style_lst.append(
        [
            {
                "chain": "A",
                "resi": 12,
            },
            {"stick": {"colorscheme": "lightgrayCarbon", "radius": 0.2}},
        ]
    )

    for loop_name, loop_resids in loop_resid_dict.items():

        if loop_name == sw1_name:
            loop_color = sw1_color
            stick_resid = 32
        elif loop_name == sw2_name:
            loop_color = sw2_color
            stick_resid = 71

        surface_lst.append(
            [
                {"opacity": opacity, "color": loop_color},
                {"chain": "A", "resi": loop_resids, "hetflag": False},
            ]
        )

        style_lst.append(
            [
                {
                    "chain": "A",
                    "resi": [loop_resids],
                },
                {
                    "cartoon": {
                        "style": "oval",
                        "color": loop_color,
                        "thickness": 0.2,
                    }
                },
            ]
        )

        style_lst.append(
            [
                {
                    "chain": "A",
                    "resi": [stick_resid],
                    "elem": "C",
                },
                {"stick": {"color": loop_color, "radius": 0.2}},
            ]
        )

        style_lst.append(
            [
                {"chain": "A", "resi": [stick_resid], "elem": ["O", "N", "H"]},
                {"stick": {"colorscheme": "whiteCarbon", "radius": 0.2}},
            ]
        )

    left_col, right_col = st.columns(2)

    with left_col:
        show_st_structure(
            "6oim",
            style_lst=style_lst,
            surface_lst=surface_lst,
            cartoon_style="oval",
            spin_on=True,
            zoom_dict={"chain": "A"},
            zoom=1.2,
            width=400,
            height=300,
        )

    right_col.markdown("# Rascore")
    right_col.markdown("### A tool for analyzing RAS protein structures")
    right_col.markdown("**Created by Mitchell Parker and Roland Dunbrack**")
    right_col.markdown("**Fox Chase Cancer Center**")

    database_link_dict = {
        "bioRxiv Paper": "https://www.biorxiv.org/content/10.1101/2022.02.02.478568v1",
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
        *Rascore* is a tool for analyzing structures of the RAS protein:
        KRAS, NRAS, and HRAS. The *Rascore* 
        database presents a continually updated analysis of all
        RAS structures in the Protein Data Bank (PDB) with their catalytic SW1 and 
        SW2 loops conformationally classified and their molecular 
        contents annotated (e.g., mutation status, nucleotide state, 
        bound protein, inhibitor site, and others). 

        Details of our work are 
        provided on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.02.478568v1) in the manuscript, 
        **An expanded classification of active, inactive and druggable RAS 
        conformations.** We hope that researchers will use 
        *Rascore* to gain novel insights into RAS biology and 
        drug discovery. 
        """
    )

    st.markdown(
        """
        ### Usage

        To the left, is a main menu for navigating to 
        each page in the *Rascore* database:

        - **Home Page:** We are here!
        - **Database Overview:** Overview of the *Rascore* database, molecular annotations, 
        and conformational classification.
        - **Explore Conformations:** Explore RAS conformations by nucleotide state.
        - **Search PDB:** Search for individual PDB entries containing RAS structures.
        - **Query Database:** Query the *Rascore* database by RAS conformations and molecular annotations.
        - **Compare Inhibitors:** Compare inhibitor-bound RAS structures by individual PDB entries or conformations.
        - **Classify Structures:** Conformationally classify your RAS structures.
        """
    )
    st.markdown("---")

    st.markdown(
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

    st.markdown("---")

    st.markdown(
        """
        ### Funding

        - NIH F30 GM142263 (to M.P.)
        - NIH R35 GM122517 (to R.D.)
         """
    )

    st.markdown(
        """
        ### License
        Apache License 2.0
        """
    )

    write_st_end()
