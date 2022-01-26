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


def home_page():

    st.sidebar.markdown("## External Links")

    st.sidebar.markdown(
        f"""
        - [BioArxiv Aricle](LINK)
        - [GitHub Page](https://github.com/mitch-parker/rascore)
        - [Dunbrack Lab](https://dunbrack.fccc.edu/retro/)
        - [RCSB Protein Data Bank](https://www.rcsb.org)
        """
    )

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
        provided on [BioArxiv](LINK) in the manuscript: **An expanded 
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
        """
        ### Authors
        Please feel free to contact us with any issues, comments, or questions.

        ##### Mitchell Parker

        - Email: <mitch.isaac.parker@gmail.com>
        - GitHub: https://github.com/mitch-parker

        ##### Roland Dunbrack

        - Email: <roland.dunbrack@gmail.com>
        - GitHub: https://github.com/DunbrackLab
        """
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

        Copyright (c) 2022 Mitchell Isaac Parker
        """
    )
