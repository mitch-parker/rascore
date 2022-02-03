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

from ..functions.gui import load_st_table


def conformation_page():

    st.markdown("# Browse Conformations")

    st.markdown("---")

    left_col, right_col = st.columns(2)

    left_col.markdown("#### Overview")
    left_col.markdown(
        """
    - We analyzed all 699 available human KRAS, NRAS, and HRAS structures in the Protein Data Bank (PDB) to 
    define a more comprehensive classification of active, inactive, and druggable RAS conformations. 
    - We first annotated the molecular contents of each RAS structure, including their mutation status, 
    nucleotide state and bound protein (e.g., effector, GAP, GEF) or inhibitor site (e.g., SP12, SP2). 
    - Second, we conformationally clustered these structures based on the arrangement of their catalytic 
    switch 1 (SW1) and switch 2 (SW2) loops to create a biologically and therapeutically informed map 
    of the RAS conformational landscape. 
    """
    )

    right_col.markdown("#### Key Terms")
    right_col.markdown(
        """

    - Bound Protein
        - **Effector:** RAS-binding domain (RBD) or RAS-associating domain (RA) containing protein.
        - **GEF.CDC25:** Guanine exchange factor catalaytic domain
            - Removes GDP allowing subsequent GTP rebinding
        - **GEF.REM:** Guanine exchange factor allosteric domain
            - Accelerates GDP release at the CDC25 domain of SOS1
        - **GAP:** GTPase activating protein
            - Catalyze the otherwise slow intrinsic rate of GTP hydrolysis to GDP
        - **Binder:** Designed protein inhibitor (e.g., Affimer, DARPin)
        - **Nanodisc:** Synthetic membrane

    - Inhibitor Site
        - **SP12:** SW1/SW2 pocket
        - **SP2:** SW2 pocket
    """
    )

    st.markdown("#### Conformational Clustering")

    st.markdown(
        """
    -	We separately clustered the arrangements of SW1 (residues 25-40) and SW2 (residues 56-76)
    based on their backbone dihedral angle values: φ (phi), ψ (psi), and ω (omega).
    - In our analysis, we used the Density-Based Spatial Clustering of Applications with Noise (DBSCAN) 
    algorithm with a distance metric that locates the maximum backbone dihedral difference upon pairwise 
    comparison of loop residues.
    - We first separated RAS structures by their nucleotide state (“0P” for nucleotide-free, “2P” for GDP-bound, 
    and “3P” for GTP or GTP analog-bound) and subsequently clustered the conformations of SW1 and SW2 for each 
    nucleotide state using DBSCAN. 
    - For clarity and brevity in our classification, we named each SW1 and SW2 conformational cluster by its loop 
    name and nucleotide state and then added further designations as needed.
    """
    )

    df = load_st_table(__file__)
