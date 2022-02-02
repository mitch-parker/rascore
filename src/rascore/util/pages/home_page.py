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
from PIL import Image

from ..functions.gui import (
    write_st_end,
    create_st_button,
)
from ..functions.path import (
    get_neighbor_path,
    get_file_path,
    pages_str,
    data_str,
)


def home_page():

    img = Image.open(
        get_file_path(
            "rascore_logo.png",
            dir_path=get_neighbor_path(__file__, pages_str, data_str),
        )
    )
    st.image(img)

    st.sidebar.markdown("## Database Links")

    database_link_dict = {
        "bioArxiv Paper": "LINK",
        "GitHub Page": "https://github.com/mitch-parker/rascore",
        "RCSB Protein Data Bank": "https://www.rcsb.org",
    }

    for link_text, link_url in database_link_dict.items():
        create_st_button(link_text, link_url, st_col=st.sidebar)

    st.sidebar.markdown("## Software Links")

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
        catalytic switch 1 (SW1) and switch 2 (SW2) loops. 
        In addition, *rascore* can be used to build and query 
        an updatable database of all available RAS structures 
        from the Protein Data Bank with their SW1 and SW2 loops 
        conformationally classified and their molecular contents 
        annotated (e.g., mutation status, nucleotide state, bound 
        proteins, small molecule inhibitors, etc.). 

        Details of our RAS conformational classification are 
        provided on [bioArxiv](LINK) in the manuscript: **An expanded 
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
        each page in the *rascore* GUI application:

        - Home Page
        - Search PDB
        - Query Database
        - Classify Structures
        """
    )
    st.markdown("---")

    st.markdown(
        f"""
        ### Authors
        Please feel free to contact us with any issues, comments, or questions.

        ##### Mitchell Parker [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40Mitch_P)](https://twitter.com/Mitch_P)

        - Email: <mitch.isaac.parker@gmail.com>
        - GitHub: https://github.com/mitch-parker

        ##### Roland Dunbrack [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social&label=Follow%20%40RolandDunbrack)](https://twitter.com/RolandDunbrack)

        - Email: <roland.dunbrack@gmail.com>
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
