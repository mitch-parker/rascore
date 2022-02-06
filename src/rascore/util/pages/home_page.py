# -*- coding: utf-8 -*-
"""
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
"""

import streamlit as st

from ..constants.pharm import sp2_color
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
            {"stick": {"color": sp2_color, "radius": 0.2}},
        ]
    )
    style_lst.append(
        [
            {
                "chain": "A",
                "resn": "MOV",
                "elem": ["N", "O", "H"],
            },
            {"stick": {"colorscheme": "Carbon", "radius": 0.2}},
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
                {"stick": {"colorscheme": "Carbon", "radius": 0.2}},
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

    right_col.markdown("# rascore")
    right_col.markdown("### A tool for analyzing the conformations of RAS structures")
    right_col.markdown("**Created by Mitchell Parker and Roland Dunbrack**")
    right_col.markdown("**Fox Chase Cancer Center**")

    st.sidebar.markdown("## Database-Related Links")

    database_link_dict = {
        "bioRxiv Paper": "LINK",
        "GitHub Page": "https://github.com/mitch-parker/rascore",
        "RCSB Protein Data Bank": "https://www.rcsb.org",
    }

    for link_text, link_url in database_link_dict.items():
        create_st_button(link_text, link_url, st_col=st.sidebar)

    st.sidebar.markdown("## Software-Related Links")

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

    left_link_col, middle_link_col, right_link_col = st.sidebar.columns(3)

    i = 0
    link_col_dict = {0: left_link_col, 1: middle_link_col, 2: right_link_col}
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
        *rascore* is a tool for analyzing RAS structures 
        (KRAS, NRAS, and HRAS) by the conformations of their 
        catalytic switch 1 (SW1) and switch 2 (SW2) loops. The *rascore*  
        database presents an updated analysis of all available 
        RAS structures in the Protein Data Bank (PDB) with their SW1 and 
        SW2 loops conformationally classified and their molecular 
        contents annotated (e.g., mutation status, nucleotide state, 
        bound protein, inhibitor site, etc.). 

        Details of our RAS conformational classification are 
        provided on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.02.478568v1) in our paper: **An expanded 
        classification of active, inactive and druggable RAS 
        conformations.** We hope that researchers will use 
        *rascore* to gain novel insights into RAS biology and 
        drug discovery. 
        """
    )

    st.markdown(
        """
        ### Usage

        To the left, there is a main menu for navigating to 
        each page in the *rascore* application:

        - **Home Page:** We are here!
        - **Database Overview:** Overview of the *rascore* database, molecular annotations, 
        and conformational classification.
        - **Explore Conformations:** Explore RAS conformations by nucleotide state.
        - **Search Database:** Search for individual PDB entries.
        - **Query Database:** Query the *rascore* database by RAS conformations and molecular annotations.
        - **Classify Structures:** Conformationally classify your inputted RAS structures.
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
        MIT License
        """
    )

    write_st_end()
