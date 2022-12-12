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

    st.markdown("## Database Overview")

    st.markdown("---")

    st.sidebar.markdown(
        "**Note.** See our [*Cancer Research*](https://aacrjournals.org/cancerres/article/doi/10.1158/0008-5472.CAN-22-0804/696349/Delineating-The-RAS-Conformational-LandscapeThe) paper for additional details."
    )

    text_col, leg_1_col, leg_2_col = st.columns((2, 0.75, 0.75))

    text_col.markdown("#### Introduction to the Rascore Database")
    text_col.markdown(
        f"""
    - For many human cancers and tumor-associated diseases, mutations in the RAS isoforms (KRAS, NRAS, and HRAS) are the most common oncogenic alterations, making these proteins high-priority therapeutic targets. 
    - Effectively targeting the RAS isoforms requires an exact understanding of their active, inactive, and druggable conformations.  
    - In consequence, we analyzed over 700 available human KRAS, NRAS, and HRAS structures in the Protein Data Bank (PDB) to create a comprehensive map of the RAS conformational landscape. 
      - First, we annotated the molecular contents of each RAS structure, including their mutation status, nucleotide state, and bound protein or inhibitor site (*see right for key terms*). 
      - Second, we conformationally classified all available RAS structures based on the configurations of their catalytic switch 1 (SW1) and switch 2 (SW2) loops, identifying three SW1 and nine SW2 conformations.
    - The *Rascore* database presents a continually updated dataset of RAS structures in the PDB that are conformationally classified and annotated for their molecular contents (*now {len(df)} structures*).
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

    st.markdown("#### Our Conformational Classification Algorithm")

    st.markdown(
        """
    - First, we broadly classified RAS structures based on the spatial positions of residue Y32 relative to the active site (“Y32in” or “Y32out”) in SW1 and residue Y71 relative to the hydrophobic core (“Y71in” or “Y71out”) in SW2. 
    - Second, we separately clustered the configurations of SW1 (residues 25-40) and SW2 (residues 56-76) based on their backbone dihedral angle values: φ (phi), ψ (psi), and ω (omega).
    - In clustering RAS SW1 and SW2 loops, we used the Density-Based Spatial Clustering of Applications with Noise (DBSCAN) algorithm with a distance metric that locates the maximum backbone dihedral angle difference upon pairwise comparison of loop residues (previously implemented by our group for other proteins, such as [kinases](http://dunbrack.fccc.edu/kincore/home) and [antibodies](http://dunbrack2.fccc.edu/PyIgClassify/)).
      - DBSCAN clusters points with sufficient numbers of near neighbors and classifies the remainder as outliers. 
      - We first separated RAS structures by nucleotide state (0P, 2P, and 3P) and spatial class (Y32in/out for SW1 and Y71in/out for SW2) and, within each group, subsequently clustered the conformations of completely modeled SW1 and SW2 loops with well-defined electron density. 
      - We then assigned a small number of poorly or incompletely modeled loops to the clusters obtained from DBSCAN by using a nearest neighbors (NN) approach that we developed.
      - In the Rascore database, we use our NN approach to conformationally classify user uploaded structures and additional RAS structures deposited to the PDB.
    """
    )

    st.markdown("---")

    st.markdown("#### Summary of the Current RAS Conformational Landscape")

    st.markdown("""
    - For clarity and brevity in our RAS conformational classification, we named each SW1 and SW2 conformation by its *spatial class*, *nucleotide state*, and a *conformational label* (written in all-capital letters). 
    - The SW1 conformations are labeled Y32in.3P-ON (GTP-bound state 2), Y32out.2P-OFF (GDP-bound), and Y32out.0P-GEF (nucleotide-free) (*panel* **A**).
    - There was no structurally uniform cluster within Y32out.3P structures that could be called the GTP-bound state 1. 
    - The only nucleotide-free SW2 conformation was Y32out.0P-GEF (*panel* **B**).
    - The GTP-bound SW2 conformations included Y71in.3P-R (R state) and Y71out.3P-T (T state), and two previously unclassified druggable conformations associated with inhibitors at the SP12 site, which we named Y71in.3P-SP12-A and Y71.3P-SP12-B (*panel* **C**). 
    - The four GDP-bound SW2 conformations are all named for their predominant binding partners, which consist of SP2 and SP12 inhibitors and protein binders: Y71out.2P-SP2-A, Y71out.2P-SP2-B, Y71in.2P-SP12, and Y71out.2P-BINDER (*panel* **D**).
    - In *panels* **E** and **F**, our original SW1 and SW2 conformational clustering (with NN assignments added) are displayed as Ramachandran maps per residue of each cluster.
    """)

    img = Image.open(
        get_file_path(
            "rascore_figure.png",
            dir_path=get_neighbor_path(__file__, pages_str, data_str),
        )
    )

    st.image(img, output_format="PNG")


    write_st_end()
