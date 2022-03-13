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

import pandas as pd
import streamlit as st

from ..constants.conf import sw1_color, sw2_color, sw1_name, sw2_name, loop_resid_dict
from ..functions.table import extract_int, lst_col, str_to_dict
from ..functions.lst import str_to_lst, lst_nums, res_to_lst
from ..functions.gui import (
    load_st_table,
    show_st_table,
    mask_st_table,
    show_st_structure,
    write_st_end,
    get_html_text,
    ribbon_name,
    trace_name,
    standard_name,
    aa_name,
)
from ..functions.col import (
    rename_col_dict,
    pdb_code_col,
    chainid_col,
    bio_lig_col,
    ion_lig_col,
    pharm_lig_col,
    chem_lig_col,
    mod_lig_col,
    mem_lig_col,
    gene_class_col,
    method_col,
    resolution_col,
    r_factor_col,
    space_col,
    mut_status_col,
    nuc_class_col,
    prot_class_col,
    match_class_col,
    pocket_class_col,
    interf_class_col,
    bound_prot_col,
    bound_prot_swiss_id_col,
    bound_prot_chainid_col,
    bound_lig_cont_col,
    bound_prot_cont_col,
    sw1_col,
    sw2_col,
    y32_col,
    y71_col,
    date_col,
)

sw1_resid_lst = res_to_lst(loop_resid_dict[sw1_name])
sw2_resid_lst = res_to_lst(loop_resid_dict[sw2_name])

def pdb_page():

    st.markdown("# Search PDB")

    st.markdown("---")

    df = load_st_table(__file__)

    st.sidebar.markdown("### PDB Selection")

    pdb_code = st.sidebar.selectbox(
        "Entry", [x.upper() for x in lst_col(df, pdb_code_col, unique=True)]
    )
    pdb_df = mask_st_table(df, {pdb_code_col: pdb_code.lower()})

    chainid = st.sidebar.selectbox("Chain", lst_col(pdb_df, chainid_col))
    chainid_df = mask_st_table(pdb_df, {chainid_col: chainid})

    gene_class = chainid_df.at[0, gene_class_col]
    mut_status = chainid_df.at[0, mut_status_col]

    st.markdown(f"#### PDB: [{pdb_code}](https://www.rcsb.org/structure/{pdb_code}) (Chain {chainid}) - {gene_class}({mut_status})")

    left_col, right_col = st.columns(2)

    left_col.markdown("##### General Information")
    for col in [method_col, resolution_col, r_factor_col, space_col, date_col]:
        left_col.markdown(f"**{rename_col_dict[col]}:** {chainid_df.at[0,col]}")

    right_col.markdown("##### Molecular Annotations")

    annot_df = pd.DataFrame()
    for i, col in enumerate(
        [
            nuc_class_col,
            prot_class_col,
            pocket_class_col,
            match_class_col,
            interf_class_col,
        ]
    ):
        annot_df.at[i, "Molecular Content"] = rename_col_dict[col]
        annot_df.at[i, "Annotation"] = chainid_df.at[0, col]

    show_st_table(annot_df, st_col=right_col)

    if chainid_df.at[0,bound_prot_col] != "None":
        for col in [bound_prot_col, bound_prot_swiss_id_col]:
            left_col.markdown(f"**{rename_col_dict[col]}:** {chainid_df.at[0,col]}")

    st_col_lst = st.columns(4)

    for i, col in enumerate([sw1_col, sw2_col, y32_col, y71_col]):

        if col in [sw1_col, y32_col]:
            color = sw1_color
        elif col in [sw2_col, y71_col]:
            color = sw2_color

        st_col_lst[i].markdown(f"##### {rename_col_dict[col]}")
        st_col_lst[i].markdown(get_html_text({chainid_df.at[0, col]:color}, 
        font_size="large"),unsafe_allow_html=True)

    st.markdown("---")

    left_check_col, middle_check_col, right_check_col = st.columns(3)
    left_view_col, right_view_col = st.columns([0.4,0.6])

    left_view_col.markdown("##### Viewer Settings")

    stick_resids = [32, 71]
    zoom_resids = None
    cartoon_trans = 1.0
    surface_trans = 0.0
    zoom = 1.5

    left_check_col.markdown("##### Bound Ligands")

    lig_check_dict = dict()
    for col in [
        bio_lig_col,
        ion_lig_col,
        pharm_lig_col,
        chem_lig_col,
        mod_lig_col,
        mem_lig_col,
    ]:
        lig_lst = str_to_lst(chainid_df.at[0, col])
        if "None" not in lig_lst:
            for lig in lig_lst:
                lig_check_dict[lig] = left_check_col.checkbox(f"{rename_col_dict[col]}: {lig}")

    if len(lig_check_dict.keys()) == 0:
        left_check_col.write("No bound ligands.")

    pharm_lig_lst = str_to_lst(chainid_df[pharm_lig_col].iloc[0]) 

    pharm_on = False
    if "None" not in pharm_lig_lst:
        if left_view_col.checkbox("Display Inhibitor Site"):

            pharm_on = True

            pharm_lig = left_view_col.selectbox("Select Inhibitor Site", pharm_lig_lst)
            pharm_cont_dict = str_to_dict(chainid_df[bound_lig_cont_col].iloc[0], return_int=True)

            stick_resids = pharm_cont_dict[pharm_lig]
            zoom_resids = {"resn": pharm_lig}
            cartoon_trans = 0.5
            surface_trans = 0.7
            zoom = 1.0

    middle_check_col.markdown("##### Bound Proteins")

    prot_check_dict = dict()
    for prot_lst in str_to_lst(chainid_df.at[0, bound_prot_chainid_col]):
        if "None" not in prot_lst:
            for prot in prot_lst:
                prot_check_dict[prot] = middle_check_col.checkbox(f"Chain {prot}")

    if len(prot_check_dict.keys()) == 0:
        middle_check_col.write("No bound proteins.")


    bound_prot_lst = str_to_lst(chainid_df[bound_prot_chainid_col].iloc[0]) 

    if "None" not in bound_prot_lst:
        if left_view_col.checkbox("Display Bound Protein Site"):

            if pharm_on:
                left_view_col.warning("Cannot display bound protein and inhibitor sites simultaneously.")
            else:
                bound_prot = left_view_col.selectbox("Select Inhibitor Site", bound_prot_lst)
                prot_cont_dict = str_to_dict(chainid_df[bound_prot_cont_col].iloc[0], return_int=True)

                stick_resids = prot_cont_dict[bound_prot]
                cartoon_trans = 0.5
                surface_trans = 0.0
                zoom = 1.0

    right_check_col.markdown("##### Mutation Sites")

    mut_check_dict = dict()
    for mut in str_to_lst(chainid_df.at[0, mut_status_col]):
        if mut != "WT":
            mut_check_dict[extract_int(mut)] = right_check_col.checkbox(mut)

    if len(mut_check_dict.keys()) == 0:
        right_check_col.write("Not mutated.")

    stick_resids = left_view_col.multiselect("Displayed Residues", lst_nums(1, 189), default=stick_resids)

    label_resids = left_view_col.checkbox("Label Residues", value=True)

    scheme_dict = {standard_name: False, aa_name: True}
    aa_scheme = scheme_dict[left_view_col.radio("Color Scheme", [standard_name, aa_name])]

    style_dict = {ribbon_name: "oval", trace_name: "trace"}
    cartoon_style = style_dict[
        left_view_col.radio("Cartoon Style", [ribbon_name, trace_name])
    ]

    cartoon_trans = left_view_col.slider(
        "Cartoon Transparency", min_value=0.0, max_value=1.0, value=cartoon_trans
    )

    surface_trans = left_view_col.slider(
        "Surface Transparency", min_value=0.0, max_value=1.0, value=surface_trans
    )

    spin_on = left_view_col.checkbox("Rotate Structure")

    all_chains = False
    if len(pdb_df) > 1:
        all_chains = left_view_col.checkbox("Show All Chains")


    show_st_structure(chainid_df,
                        mut_resids=list(mut_check_dict.keys()), stick_resids=stick_resids,
                        label_muts=mut_check_dict, label_resids=label_resids, 
                        label_ligs=lig_check_dict, label_prots=prot_check_dict, 
                        cartoon_style=cartoon_style,
                        cartoon_trans=cartoon_trans, 
                        surface_trans=surface_trans,
                        mut_trans=surface_trans,
                        aa_scheme=aa_scheme,
                        spin_on=spin_on,
                        all_chains=all_chains,
                        zoom_resids=zoom_resids,
                        zoom=zoom,
                        width=500,
                        height=500,
                        st_col=right_view_col)

    write_st_end()