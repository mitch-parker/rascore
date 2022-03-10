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
import numpy as np

from ..functions.chem import draw_lig_plot, get_lig_mcs, get_lig_smiles, get_lig_simi, is_lig_match
from ..functions.table import mask_equal, lst_col, mask_unequal, str_to_dict
from ..functions.col import (pocket_class_col, pdb_id_col, pdb_code_col, 
bound_lig_cont_col, chainid_col, pharm_lig_smiles_col, pharm_class_col,
pocket_lig_col, pharm_lig_col, pocket_score_col, match_class_col,
pocket_volume_col, pharm_lig_col, rename_col_dict)
from ..functions.file import pocket_table_file
from ..constants.pharm import sp2_name, sp12_name, bound_name, pharm_color_dict
from ..constants.conf import loop_resid_dict, loop_color_dict, sw1_name, sw2_name, y32_name, y71_name
from ..functions.gui import load_st_table, write_st_end, show_st_structure, show_st_table, get_html_text
from ..functions.lst import calc_simpson, res_to_lst,lst_unique, str_to_lst

import streamlit as st

sp12_bound_name = f"{sp12_name}-{bound_name}"
sp2_bound_name = f"{sp2_name}-{bound_name}"

def inhibitor_page():
    
    st.markdown("# Compare Inhibitors")

    st.markdown("---")

    df = load_st_table(__file__)
    pocket_df = load_st_table(__file__, file_name=pocket_table_file)

    pocket_df = mask_unequal(pocket_df, pocket_lig_col, "STP")

    df = mask_equal(df, pocket_class_col, [sp12_bound_name, sp2_bound_name])
    df = mask_equal(df, pdb_id_col, lst_col(pocket_df, pdb_id_col))

    left_query_col, right_query_col = st.columns(2)

    pocket_class = left_query_col.radio(rename_col_dict[pocket_class_col], [sp12_bound_name, sp2_bound_name])
    pharm_class = pocket_class.split('-')[0]

    pharm_df = mask_equal(df, pharm_class_col, pharm_class)

    for col in [match_class_col]:

        mask_lst = right_query_col.multiselect(rename_col_dict[col], 
                                    lst_col(pharm_df, col, unique=True))

        if len(mask_lst) > 0:
            pharm_df = mask_equal(df, col, mask_lst)


    with st.expander("Additional Site Query", expanded=False):

        cont_lst = st.multiselect("Residue Contacts", res_to_lst("1-189"))

        if len(cont_lst) > 0:

            cont_simi = st.slider("Residue Contact Similarity (Simpson)",  min_value=0.0, max_value=1.0, value=0.4)

            st.write("Querying for Sites")
            site_bar = st.progress(0)
            site_index_lst = list()
            for i, index in enumerate(list(pharm_df.index.values)):
                site_cont_dict = str_to_dict(pharm_df.at[index, bound_lig_cont_col], return_int=True)
                site_cont_lst = site_cont_dict[pharm_df.at[index, pharm_lig_col]]
                if calc_simpson(site_cont_lst, cont_lst) >= cont_simi:
                    site_index_lst.append(index)

                site_bar.progress((i + 1)/len(pharm_df))

            pharm_df = pharm_df.loc[site_index_lst, :]

    with st.expander("Additional Chemistry Query", expanded=False):

        smiles_str = st.text_input("SMILES Strings (Comma-Separated)")

        st.markdown("**Tip.** Use [PubChem Draw Tool](https://pubchem.ncbi.nlm.nih.gov/#draw=true) to design and edit SMILES strings.")

        if len(smiles_str) > 0:

            smiles_lst = str_to_lst(smiles_str)

            exact_name = 'Exact'
            simi_name = 'Similarity'

            match_type = st.radio("Substructure Match", [exact_name, simi_name])

            if match_type == exact_name:
                min_match = st.number_input("Minimum Matches", min_value=1, max_value=len(smiles_lst), value=1)

            elif match_type == simi_name:
                mean_name = "Mean"
                max_name = "Max"
                min_name = "Min"

                min_simi = st.slider("2D Fingerprint Similarity (Tonimoto)", min_value=0.0, max_value=1.0, value=0.1)

                if len(smiles_lst) > 1:
                    simi_method = st.radio("Similarity Method", [mean_name, max_name, min_name])

            st.write("Querying for Chemistries")
            chem_bar = st.progress(0.0)
            chem_index_lst = list()
            for i, index in enumerate(list(pharm_df.index.values)):
                chem_smiles_dict = str_to_dict(pharm_df.at[index, pharm_lig_smiles_col], return_str=True)
                chem_smiles= chem_smiles_dict[pharm_df.at[index, pharm_lig_col]][0]
                if match_type == exact_name:
                    if is_lig_match(lig=chem_smiles,matches=smiles_lst) >= min_match:
                        chem_index_lst.append(index)
                elif match_type == simi_name:
                    lig_simi_lst = [get_lig_simi(chem_smiles, x) for x in smiles_lst]
                    if len(smiles_lst) == 1:
                        lig_simi = lig_simi_lst[0]
                    else:
                        if simi_method == mean_name:
                            lig_simi = np.mean(lig_simi_lst)
                        elif simi_method == max_name:
                            lig_simi = np.max(lig_simi_lst)
                        elif simi_method == min_name:
                            lig_simi = np.min(lig_simi_lst)

                    if lig_simi >= min_simi:
                        chem_index_lst.append(index)

                chem_bar.progress((i + 1)/len(pharm_df))

            pharm_df = pharm_df.loc[chem_index_lst, :]

    for col in [pharm_lig_col]:

        mask_lst = st.multiselect(rename_col_dict[col], 
                                    lst_col(pharm_df, col, unique=True))

        if len(mask_lst) > 0:
            pharm_df = mask_equal(df, col, mask_lst)

    st.info(f"Total of {len(pharm_df)} inhibitor-bound structures based on provided query.")

    if len(pharm_df) < 2:
        st.warning("Insufficient Number of Structures Based On Selected Query")
    else:
        with st.expander(f"{pocket_class}: One-to-One Comparison", expanded=True):

            st.markdown(f"#### {pocket_class}: One-to-One Comparison")

            left_pdb_col, right_pdb_col = st.columns(2)

            pdb_id_lst = [x.upper() for x in lst_col(pharm_df, pdb_id_col)]
            
            pdb_id_1 = left_pdb_col.selectbox("PDB ID (Left)", pdb_id_lst)

            pdb_id_2 = right_pdb_col.selectbox("PDB ID (Right)", [x for x in pdb_id_lst if x != pdb_id_1])

            pdb_code_1 = pdb_id_1[:4].lower()
            chainid_1 = pdb_id_1[4:5]

            pdb_code_2 = pdb_id_2[:4].lower()
            chainid_2 = pdb_id_2[4:5]

            pdb_id_1 = f"{pdb_code_1}{chainid_1}"
            pdb_id_2 = f"{pdb_code_2}{chainid_2}"

            pharm_1_df = mask_equal(pocket_df, pdb_id_col, pdb_id_1)
            pharm_2_df = mask_equal(pocket_df, pdb_id_col, pdb_id_2)

            pocket_1_df = mask_equal(pocket_df, pdb_id_col, pdb_id_1)
            pocket_2_df = mask_equal(pocket_df, pdb_id_col, pdb_id_2)

            view_dict = {pdb_id_1: [left_pdb_col, pharm_1_df, pocket_1_df], pdb_id_2: [right_pdb_col, pharm_2_df, pocket_2_df]}

            st.sidebar.markdown("## Viewer Settings")

            show_cont = st.sidebar.checkbox("Label Contacts",value=True)

            style_dict = {"Ribbon": "oval", "Trace": "trace"}

            scheme_dict = {"Standard": False, "Amino Acid": True}
            amino_scheme = scheme_dict[
                st.sidebar.radio("Residue Colors", ["Standard", "Amino Acid"])
            ]

            cartoon_style = style_dict[
                st.sidebar.radio("Cartoon Style", ["Ribbon", "Trace"])
            ]

            cartoon_trans = st.sidebar.slider(
                "Cartoon Transparency", min_value=0.0, max_value=1.0, value=0.5
            )

            surf_trans = st.sidebar.slider(
                "Surface Transparency", min_value=0.0, max_value=1.0, value=0.0
            )

            for pdb_id in list(view_dict.keys()):
                view_col = view_dict[pdb_id][0]
                pharm_pdb_df = view_dict[pdb_id][1]
                pocket_pdb_df = view_dict[pdb_id][1]

                pdb_code = pharm_pdb_df[pdb_code_col].iloc[0]
                chainid = pharm_pdb_df[chainid_col].iloc[0]
                lig = pharm_pdb_df[pharm_lig_col].iloc[0]
                match_class = pharm_pdb_df[match_class_col].iloc[0]
                smiles_dict = str_to_dict(pharm_pdb_df[pharm_lig_smiles_col].iloc[0], return_str=True)
                cont_dict = str_to_dict(pharm_pdb_df[bound_lig_cont_col].iloc[0], return_int=True)

                smiles = smiles_dict[lig][0]

                view_col.markdown(f"**{rename_col_dict[pharm_lig_col]}:** [{lig}](https://www.rcsb.org/ligand/{lig}) ({match_class})")
                view_col.markdown(f"**SMILES:** {smiles}")
                
                for col in [pocket_score_col, pocket_volume_col]:
                    view_col.markdown(f"**{rename_col_dict[col]}:** {round(float(pocket_pdb_df[col].max()),2)}")

                style_lst = list()
                reslabel_lst = list()
                surface_lst = list()

                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "invert": True,
                        },
                        {
                            "cartoon": {
                                "color": "white",
                                "style": cartoon_style,
                                "thickness": 0.2,
                                "opacity": 0,
                            }
                        },
                    ]
                )

                style_lst.append(
                    [
                        {
                            "chain": chainid,
                        },
                        {
                            "cartoon": {
                                "color": "white",
                                "style": cartoon_style,
                                "thickness": 0.2,
                                "opacity": cartoon_trans,
                            }
                        },
                    ]
                )  

                style_lst.append(
                        [
                            {
                                "chain": chainid,
                                "resn": [lig],
                                "elem": "C",
                            },
                            {"stick": {"color": pharm_color_dict[pharm_class], "radius": 0.2}},
                        ]
                    )

                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "resn": [lig],
                            "elem": ["N", "O", "H"],
                        },
                        {"stick": {"colorscheme": "lightgrayCarbon", "radius": 0.2}},
                    ]
                )

                if sp2_name in pocket_class:
                    style_lst.append(
                        [
                            {
                                "chain": chainid,
                                "resi": 12,
                            },
                            {"stick": {"colorscheme": "lightgrayCarbon", "radius": 0.2}},
                        ]
                    )

                surface_lst = [
                [
                    {"opacity": surf_trans, "color": "white"},
                    {"chain": chainid, "hetflag": False},
                ]
            ]

                sw1_resid_lst = res_to_lst(loop_resid_dict[sw1_name])
                sw2_resid_lst = res_to_lst(loop_resid_dict[sw2_name])

                if amino_scheme:
                    style_lst.append(
                        [
                            {"chain": chainid,
                                "resi": cont_dict[lig],
                                "elem": "C",},
                            {"stick": {"colorscheme": 'amino', "radius": 0.2}},
                        ]
                    ) 
                    style_lst.append(
                            [
                                {"chain": chainid, "resi": cont_dict[lig], "elem": ["O", "N", "H"]},
                                {"stick": {"colorscheme": "whiteCarbon", "radius": 0.2}},
                            ]
                        )

                else:
                    for resid in cont_dict[lig]:

                        resid_color = "white"
                        if int(resid) in sw1_resid_lst:
                            resid_color = loop_color_dict[sw1_name]
                        if int(resid) in sw2_resid_lst:
                            resid_color = loop_color_dict[sw2_name]


                        style_lst.append(
                        [
                            {
                                "chain": chainid,
                                "resi": [resid],
                                "elem": "C",
                            },
                            {"stick": {"color": resid_color, "radius": 0.2}},
                        ]
                    )

                        style_lst.append(
                            [
                                {"chain": chainid, "resi": [resid], "elem": ["O", "N", "H"]},
                                {"stick": {"colorscheme": "whiteCarbon", "radius": 0.2}},
                            ]
                        )

                if show_cont:

                    reslabel_lst.append(
                            [
                                {"chain": chainid,
                                "resi": cont_dict[lig]},
                                {
                                    "backgroundColor": "lightgray",
                                    "fontColor": "black",
                                    "backgroundOpacity": 0.5,
                                },
                            ]
                        )

                for loop_name, loop_resids in loop_resid_dict.items():

                    loop_color = loop_color_dict[loop_name]

                    surface_lst.append(
                    [
                        {"opacity": surf_trans, "color": loop_color},
                        {"chain": chainid, "resi": loop_resids, "hetflag": False},
                    ]
                )

                    style_lst.append(
                        [
                            {
                                "chain": chainid,
                                "resi": [loop_resids],
                            },
                            {
                                "cartoon": {
                                    "style": cartoon_style,
                                    "color": loop_color,
                                    "thickness": 0.2,
                                    "opacity": cartoon_trans,
                                }
                            },
                        ]
                    )

                    

                with view_col:
                    show_st_structure(
                        pdb_code,
                        style_lst=style_lst,
                        reslabel_lst=reslabel_lst,
                        surface_lst=surface_lst,
                        cartoon_style=cartoon_style,
                        zoom_dict={"chain": chainid, "resn": lig},
                        zoom=0.7,
                        width=450,
                        height=300,
                    )

                label_size = "medium"

                for col in [sw1_name, sw2_name, y32_name, y71_name]:

                    label_str = get_html_text({f"{rename_col_dict[col]}: ": "#31333F"}, font_weight='bold', font_size=label_size)

                    if col in [y32_name, sw1_name]:
                        label_color = loop_color_dict[sw1_name]
                    elif col in [y71_name, sw2_name]:
                        label_color = loop_color_dict[sw2_name]

                    label_str += get_html_text({pharm_pdb_df[col].iloc[0]: label_color}, font_size=label_size)

                    view_col.markdown(label_str, unsafe_allow_html=True)

            
            if amino_scheme:

                st.sidebar.markdown("---")

                st.sidebar.markdown("## Amino Acid Colors")

                aa_type_dict = {'Aromatic':['PHE', 'TYR', 'TRP'],
                'Acidic':['ASP','GLU'],
                'Basic':['LYS', 'ARG','HIS'],
                'Nonpolar':['ILE', 'VAL', 'LEU', 'MET','PRO', 'GLY', 'ALA'],
                'Polar':['ASN','GLN','SER','THR', 'CYS']}

                aa_color_dict = {'ASP':'#E60A0A','GLU':'#E60A0A','CYS':'#E6E600','MET':'#E6E600',
                                'LYS':'#145AFF','ARG':'#145AFF','SER':'#FA9600','THR':'#FA9600',
                                'PHE':'#3232AA','TYR':'#3232AA','ASN':'#00DCDC','GLN':'#00DCDC',
                                'GLY':'#C8C8C8','LEU':'#0F820F','VAL':'#0F820F','ILE':'#0F820F',
                                'ALA':'#C8C8C8','TRP':'#B45AB4','HIS':'#8282D2','PRO':'#DC9682'}


                aa_font_size = "medium"

                for aa_type, aa_lst in aa_type_dict.items():

                    aa_str = get_html_text({f"{aa_type}: ": "#31333F"}, font_weight='bold', font_size=aa_font_size)

                    for i, aa_name in enumerate(aa_lst): 

                        aa_str += get_html_text({aa_name:aa_color_dict[aa_name]}, font_size=aa_font_size)

                        if i != len(aa_lst) - 1:
                            aa_str += get_html_text({", ":"#31333F"}, font_size=aa_font_size)
        
                    st.sidebar.markdown(aa_str, unsafe_allow_html=True)


            st.markdown("---")

            left_chem_col, right_chem_col = st.columns(2)

            lig_1 = pharm_1_df[pharm_lig_col].iloc[0]
            lig_2 = pharm_2_df[pharm_lig_col].iloc[0]

            smiles_dict_1 = str_to_dict(pharm_1_df[pharm_lig_smiles_col].iloc[0])
            smiles_dict_2 = str_to_dict(pharm_2_df[pharm_lig_smiles_col].iloc[0])
            
            smiles_1 = smiles_dict_1[lig_1][0]
            smiles_2 = smiles_dict_2[lig_2][0]

            right_chem_col.markdown("#### Chemistry Comparison")

            mcs = get_lig_mcs([smiles_1, smiles_2])

            right_chem_col.markdown(f"**2D Fingerprint Similarity (Tanimoto):** {round(get_lig_simi(smiles_1, smiles_2),2)}")

            right_chem_col.markdown(f"**Maximum Common Substructure (MCS):** {get_lig_smiles(mcs)}")

            if not right_chem_col.checkbox(f"Show MCS",value=True):
                mcs = None
            
            right_chem_col.image(draw_lig_plot([smiles_1, smiles_2], 
            lig_labels=[f"{pdb_id_1.upper()} ({lig_1})", f"{pdb_id_2.upper()} ({lig_2})"], font_size=16,
            plot_height=4, plot_width=3, highlight_querys=mcs, color_palette=[pharm_color_dict[pharm_class]]))

            cont_dict_1 = str_to_dict(pharm_1_df[bound_lig_cont_col].iloc[0])
            cont_dict_2 = str_to_dict(pharm_2_df[bound_lig_cont_col].iloc[0])

            cont_lst_1 = cont_dict_1[lig_1]
            cont_lst_2 = cont_dict_2[lig_2]
            
            cont_df = pd.DataFrame()

            left_chem_col.markdown("#### Site Comparison")

            left_chem_col.markdown(f"**Residue Contact Similarity (Simpson):** {round(calc_simpson(cont_lst_1, cont_lst_2),2)}")

            cont_col_1 = f"{pdb_id_1.upper()} ({lig_1})"
            cont_col_2 = f"{pdb_id_2.upper()} ({lig_2})"

            for resid in sorted(lst_unique(cont_lst_1 + cont_lst_2,return_int=True)):

                cont_1 = ''
                cont_2 = ''

                if str(resid) in cont_lst_1:
                    cont_1 = "✓"
                if str(resid) in cont_lst_2:
                    cont_2 = "✓"

                cont_df.at[resid, cont_col_1] = cont_1
                cont_df.at[resid, cont_col_2] = cont_2

            cont_df = cont_df.rename_axis('Residue Contact').reset_index()

            show_st_table(cont_df, st_col=left_chem_col)

    write_st_end()