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

from ..constants.conf import (
    loop_color_dict,
)
from ..functions.col import (
    rename_col_dict,
    nuc_class_col,
    pocket_class_col,
    prot_class_col,
)
from ..constants.nuc import nuc_color_dict
from ..constants.pharm import pharm_color_dict, sp2_name, sp12_name
from ..constants.prot import (
    prot_color_dict,
    effector_name,
    gef_cdc_name,
    gef_rem_name,
    gap_name,
    binder_name,
    nano_name,
)
from ..functions.gui import (
    load_st_table,
    get_html_text,
    write_st_end,
    get_neighbor_path,
)
from ..functions.path import get_file_path, pages_str, data_str


def overview_page():

    df = load_st_table(__file__)

    st.markdown("# Database Overview")

    st.markdown("---")

    st.sidebar.markdown(
        "**Note.** See our [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.02.478568v1) paper for additional details."
    )

    text_col, leg_1_col, leg_2_col = st.columns((2, 0.75, 0.75))

    text_col.markdown("#### Building the Rascore Database")
    text_col.markdown(
        f"""
    - RAS (KRAS, NRAS, and HRAS) proteins have widespread command of cellular circuitry and are high-priority 
    drug targets in cancers and other diseases. 
    - Effectively targeting RAS proteins requires an exact understanding of their active, inactive, and 
    druggable conformations, and the structural impact of mutations. 
    - We analyzed all 699 available human KRAS, NRAS, and HRAS structures in the Protein Data Bank (PDB) to 
    define a more comprehensive mapping of the RAS conformational landscape. 
      - First, we annotated the molecular contents of each RAS structure, including their mutation status, 
    nucleotide state and bound protein or inhibitor site (*see right for key terms*). 
      - Second, we conformationally clustered all RAS structures based on the arrangement of their catalytic 
    switch 1 (SW1) and switch 2 (SW2) loops, identifying three SW1 and nine SW2 conformations (*see below for details*).
    - The *rascore* database presents a continually updated dataset of annotated and conformationally classified 
    RAS structures from the PDB (*now {len(df)} structures*).
    """
    )

    loop_term_dict = {"SW1": "Switch 1", "SW2": "Switch 2"}

    nuc_term_dict = {
        "0P": "Nucleotide-free",
        "2P": "GDP-bound",
        "3P": "GTP or GTP analog-bound",
    }

    prot_term_dict = {
        effector_name: "RAS-binding domain (RBD) or RAS-associating domain (RA) protein",
        gef_cdc_name: "Guanine exchange factor, catalaytic domain",
        gef_rem_name: "Guanine exchange factor, allosteric domain",
        gap_name: "GTPase activating protein",
        binder_name: "Designed protein inhibitor",
        nano_name: "Synthetic membrane",
    }

    pharm_term_dict = {
        sp12_name: "SW1/SW2 pocket",
        sp2_name: "SW2 pocket",
    }

    loop_str = "**Loop Names**"
    nuc_str = f"**{rename_col_dict[nuc_class_col]}**"
    pharm_str = f"**{rename_col_dict[pocket_class_col]}**"
    prot_str = f"**{rename_col_dict[prot_class_col]}**"

    leg_info_dict = {
        loop_str: loop_term_dict,
        nuc_str: nuc_term_dict,
        pharm_str: pharm_term_dict,
        prot_str: prot_term_dict,
    }
    leg_color_dict = {
        loop_str: loop_color_dict,
        nuc_str: nuc_color_dict,
        pharm_str: pharm_color_dict,
        prot_str: prot_color_dict,
    }
    leg_col_dict = {
        loop_str: leg_1_col,
        nuc_str: leg_1_col,
        pharm_str: leg_1_col,
        prot_str: leg_2_col,
    }

    for info in list(leg_info_dict.keys()):

        col = leg_col_dict[info]

        col.markdown(info)

        for term, desc in leg_info_dict[info].items():
            col.markdown(
                get_html_text(
                    {
                        term: leg_color_dict[info][term],
                        " - ": "#31333F",
                        desc: "#31333F",
                    },
                    font_size="small",
                ),
                unsafe_allow_html=True,
            )

    st.markdown("---")

    st.markdown("#### Our Conformational Clustering Algorithm")

    st.markdown(
        """
    - We separately clustered the arrangements of SW1 (residues 25-40) and SW2 (residues 56-76)
    based on their backbone dihedral angle values: φ (phi), ψ (psi), and ω (omega).
    - In our analysis, we used the Density-Based Spatial Clustering of Applications with Noise (DBSCAN) 
    algorithm with a distance metric that locates the maximum backbone dihedral difference upon pairwise 
    comparison of loop residues (previously implemented by our group for other proteins).
      - DBSCAN finds major clusters and removes outliers (labeling them as “noise”).
      - We first separated RAS structures by their nucleotide state (0P, 2P, and 3P) and subsequently 
    conformationally clustered the well modeled (complete with high electron density scores) SW1 and SW2 
    loops for each nucleotide state with DBSCAN.
      - We then assigned a small number of poorly or incompletely modeled loops to the clusters obtained from 
    DBSCAN by using a nearest neighbors (NN) approach.
      - In the rascore database, we use our NN approach to conformationally classify additional RAS structures
    deposited to the PDB.
    """
    )

    st.markdown("---")

    left_col, right_col = st.columns((2, 1.5))

    img = Image.open(
        get_file_path(
            "rascore_figure.png",
            dir_path=get_neighbor_path(__file__, pages_str, data_str),
        )
    )

    right_col.image(
        img,
    )

    left_col.markdown("#### Currently Identified SW1 and SW2 Conformations")

    left_col.markdown(
        """
    - For clarity and brevity in our classification, we named each SW1 and SW2 conformational 
    cluster by its loop name and nucleotide state and then added further designations as needed.
    - For the SW1 conformations, there was a one-to-one correspondence with the nucleotide state, 
    and we, therefore, labeled these conformations SW1.0P, SW1.2P, and SW1.3P (***panel a***).
      - These SW1 conformations can be differentiated by position of residue Y32 in SW1, 
    which was the original method for classifying these conformations
      - Previously, hydrolytically-relevant substates of GTP-bound RAS (SW1.3P) have been described based on 
    differences in hydrogen-bonding (HB) of the hydroxyl (OH) atom of Y32 with one of the γ-phosphate oxygen 
    (O1G) atoms of GTP: direct (hydrolytically incompetent), water-mediated (prefers intrinsic hydrolysis), 
    and absent (prefers GAP-mediated hydrolysis).
      - Examination of the distribution of distances between the Y32(OH) atom and the closest 
    3P(O1G) atom across RAS structures in the GTP-bound cluster (SW1.3P) revealed three peaks at 
    distances of 3, 4.5, and 7 Å; we associated these peaks with the GTP-bound substates, naming 
    them SW1.3P-Direct, SW1.3P-WaterHB, and SW1.3P-NoHB, respectively (***panels b and c***).
    - There were nine SW2 conformations in total (***panels d-f***): 
      - (*Clusters 1-2*) Previously described "R state" (SW2.3P-R) and "T state" (SW2.3P-T)
      - (*Cluster 3*) The SW2 conformation found in nucleotide-free structures (which we named 
    SW2.0P-GEF for its binding to GEFs)
      - (*Clusters 4-9*) Six previously unclassified druggable conformations, which we named by 
      their associated bound protein (only SW2.2P-Binder) or inhibitor site (SP12 or SP2) and, in 
        some cases, an indicator of cluster size order (A or B). 
    - The SW2 conformations are visualized with residue Y71 displayed since we demonstrated that the 
    position of this residue relates to RAS druggability. 
    """
    )

    write_st_end()
