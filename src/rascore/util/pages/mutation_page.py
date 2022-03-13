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

from ..scripts.make_facet_plot import make_facet_plot
from ..functions.color import change_hex_alpha, get_hex, gray_hex
from ..functions.chem import draw_lig_plot, get_lig_mcs, get_lig_smiles, get_lig_simi, is_lig_match
from ..functions.table import build_col_count_dict, mask_equal, lst_col, mask_unequal, str_to_dict, fix_col, make_dict
from ..functions.col import (pocket_class_col, pdb_id_col, pdb_code_col, id_col,
                            bound_lig_cont_col, chainid_col, pharm_lig_smiles_col, pharm_class_col,
                            nuc_class_col, mut_status_col, prot_class_col, interf_class_col,
                            pocket_lig_col, pharm_lig_col, pocket_score_col, match_class_col, 
                            modelid_col, gene_class_col,
                            pocket_volume_col, pharm_lig_col, rename_col_dict)
from ..functions.file import pocket_table_file, entry_table_file
from ..functions.path import rascore_str
from ..constants.pharm import sp2_name, sp12_name, bound_name, pharm_color_dict
from ..constants.conf import loop_resid_dict, loop_color_dict, sw1_name, sw2_name, y32_name, y71_name
from ..functions.gui import (load_st_table, show_st_dataframe, write_st_end,
                            show_st_structure, show_st_table, get_html_text, 
                            download_st_df, rename_st_cols, reorder_st_cols)
from ..functions.lst import calc_simpson, res_to_lst,lst_unique, str_to_lst

import streamlit as st

reverse_col_dict = make_dict(list(rename_col_dict.values()),list(rename_col_dict.keys()))

def mutation_page():
    
    st.markdown("## Analyze Mutations")

    st.markdown("---")

    df = load_st_table(__file__)

    st.sidebar.markdown("**Note.** We recommend comparing mutated structures within identical SW1 and SW2 conformations.")

    left_name = "Left"
    right_name = "Right"

    a_name = "Group A"
    b_name = "Group B"

    query_label_dict = {left_name: a_name, right_name: b_name}

    for query_name, label_name in query_label_dict.items():

        query_label_dict[query_name] = st.sidebar.text_input(label=f"{query_name} Label", value=label_name)

    left_name = query_label_dict[left_name]
    right_name = query_label_dict[right_name]

    left_query_col, right_query_col = st.columns(2)

    query_col_dict = {left_name: left_query_col, right_name: right_query_col}

    query_df_dict = dict()

    mut_status_lst = lst_col(df, mut_status_col, unique=True)

    for query_name, query_col in query_col_dict.items():

        query_df = df.copy(deep=True)

        query_label_dict[query_name] = label_name

        mut_status = query_col.selectbox(f"{rename_col_dict[mut_status_col]} ({query_name})", mut_status_lst)

        mut_status_lst.remove(mut_status)

        query_df = mask_equal(query_df, mut_status_col, mut_status)

        if query_col.checkbox(f"Display Selection Options ({query_name})", value=False):

            query_col.markdown(f"#### Conformation Selection ({query_name})")

            for col in [sw1_name, sw2_name, y32_name, y71_name]:

                mask_lst = query_col.multiselect(f"{rename_col_dict[col]} ({query_name})", lst_col(query_df, col, unique=True))

                if len(mask_lst) > 0:
                    query_df = mask_equal(query_df, col, mask_lst)

            query_col.markdown(f"#### Annotation Selection ({query_name})")

            for col in [gene_class_col, nuc_class_col, prot_class_col, pocket_class_col, match_class_col, interf_class_col]:

                mask_lst = query_col.multiselect(f"{rename_col_dict[col]} ({query_name})", lst_col(query_df, col, unique=True))

                if len(mask_lst) > 0:
                    query_df = mask_equal(query_df, col, mask_lst)

            query_df_dict[query_name] = query_df

    with st.expander("One-to-One Comparison", expanded=False):

        st.markdown("#### One-to-One")

        left_info_col, right_info_col = st.columns(2)   

        info_col_dict = {left_name: left_info_col, right_name: right_info_col}   

        pdb_name_dict = dict()
        pdb_df_dict = dict()

        for query_name, query_df in query_df_dict.items():

            pdb_id_lst = [x.upper() for x in lst_col(query_df, pdb_id_col) if x.upper() not in list(pdb_name_dict.values())]

            info_col = info_col_dict[query_name]
            
            pdb_upper = info_col.selectbox(f"PDB ID ({query_name})", pdb_id_lst)

            pdb_code = pdb_upper[:4].lower()
            chainid = pdb_upper[4:5]

            pdb_id = f"{pdb_code}{chainid}"
             
            pdb_df = mask_equal(query_df, pdb_id_col, pdb_id)
                
            pdb_name_dict[query_name] = pdb_upper
            pdb_df_dict[query_name] = pdb_df

            pdb_code = pdb_df[pdb_code_col].iloc[0]
            chainid = pdb_df[chainid_col].iloc[0]
            gene = pdb_df[gene_class_col].iloc[0]
            
            info_col.markdown(f" ##### PDB: [{pdb_code.upper()}](https://www.rcsb.org/structure/{pdb_code}) ({gene}) - Chain {chainid}")

                
            for col in [pocket_score_col, pocket_volume_col]:
                info_col.markdown(f"**{rename_col_dict[col]}:** {round(float(pdb_df[col].iloc[0]),2)}")

            st.markdown("---")

            st.markdown("##### Viewer Settings") 

            left_set_col, right_set_col = st.columns(2)  

            show_cont = left_set_col.checkbox("Label Residues",value=True)

            style_dict = {"Ribbon": "oval", "Trace": "trace"}

            scheme_dict = {"Standard": False, "Amino Acid": True}
            amino_scheme = scheme_dict[
                left_set_col.radio("Residue Color Scheme", ["Standard", "Amino Acid"])
            ]

            cartoon_style = style_dict[
                left_set_col.radio("Cartoon Style", ["Ribbon", "Trace"])
            ]

            cartoon_trans = right_set_col.slider(
                "Cartoon Transparency", min_value=0.0, max_value=1.0, value=0.5
            )

            surf_trans = right_set_col.slider(
                "Surface Transparency", min_value=0.0, max_value=1.0, value=0.7
            )

            left_view_col, right_view_col = st.columns(2)   

            view_col_dict = {left_name: left_view_col, right_name: right_view_col}  

            for query_name, query_df in query_df_dict.items():
                
                pdb_df = pdb_df_dict[query_name] 

                view_col = view_col_dict[query_name]

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
                        },
                        {"stick": {"radius": 0.2}},
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
                            {"stick": {"colorscheme": "amino", "radius": 0.2}},
                        ]
                    ) 
                    style_lst.append(
                            [
                                {"chain": chainid, "resi": cont_dict[lig]},
                                {"stick": {"radius": 0.2}},
                            ]
                        )
                    surface_lst.append(
                    [
                        {"opacity": surf_trans,"colorscheme": "amino"},
                        {"chain": chainid, "resi": cont_dict[lig], "hetflag": False},
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
                                {"chain": chainid, "resi": [resid]},
                                {"stick": {"radius": 0.2}},
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
                
                if not amino_scheme:
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
                        width=400,
                        height=300,
                    )

                label_size = "medium"

                if not amino_scheme:

                    for col in [sw1_name, sw2_name, y32_name, y71_name]:

                        label_str = get_html_text({f"{rename_col_dict[col]}: ": "#31333F"}, font_weight='bold', font_size=label_size)

                        if col in [y32_name, sw1_name]:
                            label_color = loop_color_dict[sw1_name]
                        elif col in [y71_name, sw2_name]:
                            label_color = loop_color_dict[sw2_name]

                        label_str += get_html_text({pdb_df[col].iloc[0]: label_color}, font_size=label_size)

                        view_col.markdown(label_str, unsafe_allow_html=True)

            if amino_scheme:

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
        
                    st.markdown(aa_str, unsafe_allow_html=True)

            st.markdown("---")

            left_site_col, right_chem_col = st.columns(2)

            pdb_1 = pdb_name_dict[left_name]
            pdb_2 = pdb_name_dict[right_name]

            pdb_1_df = pdb_df_dict[left_name]
            pdb_2_df = pdb_df_dict[right_name]

            lig_1 = pdb_1_df[pharm_lig_col].iloc[0]
            lig_2 = pdb_2_df[pharm_lig_col].iloc[0]

            smiles_dict_1 = str_to_dict(pdb_1_df[pharm_lig_smiles_col].iloc[0])
            smiles_dict_2 = str_to_dict(pdb_2_df[pharm_lig_smiles_col].iloc[0])
            
            smiles_1 = smiles_dict_1[lig_1][0]
            smiles_2 = smiles_dict_2[lig_2][0]

            right_chem_col.markdown("#### Chemistry Comparison")

            mcs = get_lig_mcs([smiles_1, smiles_2])

            right_chem_col.markdown(f"**2D Fingerprint Similarity:** {round(get_lig_simi(smiles_1, smiles_2),2)}")

            right_chem_col.markdown(f"**Maximum Common Substructure (MCS):** {get_lig_smiles(mcs)}")

            if not right_chem_col.checkbox("Highlight MCS",value=True):
                mcs = None
            
            right_chem_col.image(draw_lig_plot([smiles_1, smiles_2], 
            lig_labels=[f"{pdb_1} ({lig_1})", f"{pdb_2} ({lig_2})"], font_size=16, mol_pad=0,
            plot_height=4, plot_width=3, highlight_querys=mcs, color_palette=[pharm_color_dict[pharm_class]]))

            cont_dict_1 = str_to_dict(pdb_1_df[bound_lig_cont_col].iloc[0])
            cont_dict_2 = str_to_dict(pdb_2_df[bound_lig_cont_col].iloc[0])

            cont_lst_1 = cont_dict_1[lig_1]
            cont_lst_2 = cont_dict_2[lig_2]
            
            cont_df = pd.DataFrame()

            left_site_col.markdown("#### Site Comparison")

            left_site_col.markdown(f"**Residue Contact Similarity (Simpson):** {round(calc_simpson(cont_lst_1, cont_lst_2),2)}")

            cont_col_1 = f"{pdb_1} ({lig_1})"
            cont_col_2 = f"{pdb_2} ({lig_2})"

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

            show_st_table(cont_df, st_col=left_site_col)


        with st.expander("Many-to-Many Comparison", expanded=False):

            st.markdown("#### Many-to-Many Comparison")

            count_dict = build_col_count_dict(many_df, id_col)

            color_dict = {left_name: pharm_color_dict[pharm_class], 
            right_name: get_hex(change_hex_alpha(pharm_color_dict[pharm_class], 0.5)), both_name: gray_hex}

            left_plot_col, middle_plot_col, right_plot_col = st.columns(3)

            for label_name, label_color in color_dict.items():

                count = 0
                if label_name in list(count_dict.keys()):
                    count = count_dict[label_name]

                color_dict[label_name] = left_plot_col.color_picker(label=f"{label_name} (N={count})", 
                                                                value=label_color)

            many_order =  [x for x in list(color_dict.keys()) if x in lst_col(many_df, id_col, unique=True)]

            scatter_name = "Scatterplot"
            box_name = "Boxplot"

            plot_type = middle_plot_col.radio("Plot Type", [scatter_name, box_name])

            marker_size = middle_plot_col.number_input("Marker Size", min_value=1, max_value=50, value=5)
            
            if plot_type == scatter_name:
                plot_reg = middle_plot_col.checkbox("Display Regression",value=True)

                fig = make_facet_plot(many_df, 
                    x_col=pocket_volume_col,
                    y_col=pocket_score_col,
                    x_str=rename_col_dict[pocket_volume_col],
                    y_str=rename_col_dict[pocket_score_col],
                    x_round=0,
                    y_round=1,
                    hue_col=id_col,
                    hue_palette=color_dict,
                    hue_order=many_order,
                    plot_width=3,
                    plot_height=3,
                    font_size=12,
                    x_rotation=45,
                    marker_size=marker_size,
                    x_pad=25,
                    y_pad=10,
                    x_ha="right",
                    plot_reg=plot_reg,
                    log_reg=True,
                    trun_reg=True,
                    x_lim=[0, 1700],
                    x_ticks=[0, 500, 1000, 1500],
                    y_lim=[0, 1.2],
                    y_ticks=[0.0, 0.5, 1.0])
            
            elif plot_type == box_name:

                y_str = middle_plot_col.radio("Y-Axis",[rename_col_dict[pocket_score_col],rename_col_dict[pocket_volume_col]])

                y_col = reverse_col_dict[y_str]

                if y_col == pocket_score_col:
                    y_lim = [0, 1.2]
                    y_ticks = [0.0, 0.5, 1.0]
                    y_round = 1

                elif y_col == pocket_volume_col:
                    y_lim = [0, 1700]
                    y_ticks = [0, 500, 1000, 1500]
                    y_round = 0

                stat_pairs = None
                if len(many_order) > 1:
                    if middle_plot_col.checkbox("Display Welch's t-test p-value",value=True):
                        stat_pairs = list()
                        for id_1 in many_order:
                            for id_2 in many_order:
                                if id_1 != id_2:
                                    if (id_2, id_1) not in stat_pairs:
                                        stat_pairs.append((id_1, id_2))

                fig = make_facet_plot(
                        many_df,
                        x_col=id_col,
                        x_order=many_order,
                        x_palette=color_dict,
                        x_count=True,
                        x_str="",
                        y_col=y_col,
                        y_str=y_str,
                        show_legend=False,
                        plot_width=3,
                        plot_height=3,
                        plot_kind="box",
                        font_size=12,
                        marker_size=marker_size,
                        x_pad=35,
                        y_lim=y_lim,
                        y_ticks=y_ticks,
                        y_round=y_round,
                        stat_pairs=stat_pairs,
                        stat_test="t-test_welch",
                        stat_format="simple"
                
                    
                    )

            right_plot_col.pyplot(fig)

            many_df[pocket_score_col] = many_df[pocket_score_col].map(float)
            many_df[pocket_volume_col] = many_df[pocket_volume_col].map(float)

            left_stat_col, right_stat_col = st.columns(2)

            stat_col_dict = {pocket_score_col: left_stat_col, pocket_volume_col: right_stat_col}

            for stat_name, stat_col in stat_col_dict.items():

                stat_col.markdown(f"##### {rename_col_dict[stat_name]}")

                temp_df = pd.pivot_table(many_df, values=stat_name,
                                        index=id_col,
                                        aggfunc=[np.mean, min, max])

                stat_df = pd.DataFrame()
                for index in list(temp_df.index.values):
                    for col, stat in zip(list(temp_df.columns), [mean_name, min_name, max_name]):
                        stat_df.at[index, stat] = temp_df.at[index, col]

                for col in list(stat_df.columns):
                    stat_df[col] = stat_df[col].round(1)
                    stat_df[col] = stat_df[col].map(str)
                    stat_df = fix_col(stat_df, col)

                stat_df = stat_df.loc[many_order, :]

                show_st_table(stat_df.rename_axis(rename_col_dict[id_col]).reset_index(), st_col=stat_col)

            st.markdown("---")

            cont_df = pd.DataFrame()

            cont_order = list()

            for index in list(many_df.index.values):
                
                lig = many_df.at[index, pharm_lig_col]
                cont_dict = str_to_dict(many_df.at[index, bound_lig_cont_col], return_int=True)

                cont_lst = cont_dict[lig]

                label = many_df.at[index, id_col]
                label_name = f"{label} (N={count_dict[label]})"

                if label_name not in list(cont_df.columns):
                    cont_df[label_name] = 0

                for cont in cont_lst:
                    if cont in cont_order:
                        cont_df.at[cont, label_name] += 1
                    else:
                        cont_df.at[cont, label_name] = 1
                        cont_order.append(cont)

            cont_df = cont_df.fillna(0)

            for col in list(cont_df.columns):
                cont_df[col] = (cont_df[col] / count_dict[col.split(" (")[0]]) * 100
                cont_df[col] = cont_df[col].round(1)
                cont_df[col] = cont_df[col].map(str)
                cont_df = fix_col(cont_df, col)
                cont_df[col] += "%"
                cont_df[col] = cont_df[col].replace("0%","")

            cont_df = cont_df.loc[sorted(cont_order), [f"{x} (N={count_dict[x]})" for x in many_order if x in list(count_dict.keys())]]

            cont_df = cont_df.rename_axis('Residue Contact').reset_index()

            left_table_col, right_table_col = st.columns(2)

            right_table_col.markdown("##### Site Comparison")

            show_st_table(cont_df, st_col=right_table_col)

            for table_col in [y32_name, y71_name, sw1_name, sw2_name]:

                left_table_col.markdown(f"##### {rename_col_dict[table_col]}")

                loop_df = (
                    pd.pivot_table(
                        data=rename_st_cols(many_df),
                        index=[rename_col_dict[table_col]],
                        columns=rename_col_dict[id_col],
                        values=rename_col_dict[pdb_id_col],
                        aggfunc="nunique",
                        margins=True,
                    )
                    .fillna("")
                )

                for col in list(loop_df.columns):
                    loop_df[col] = loop_df[col].map(str)
                    loop_df = fix_col(loop_df, col)

                loop_df = reorder_st_cols(loop_df, table_col, id_col)

                loop_df = loop_df.reset_index()

                show_st_table(loop_df, st_col=left_table_col)

        #     right_table_col.markdown("##### Chemistry Comparison")

        # index_lst = list(many_df.index.values)

        # chem_df = pd.DataFrame()
        # compare_bar = st.progress(0)
        # for i, index_1 in enumerate(index_lst):
        #     smiles_1_dict = str_to_dict(many_df.at[index_1, pharm_lig_smiles_col], return_str=True)
        #     smiles_1 = smiles_1_dict[many_df.at[index_1, pharm_lig_col]][0]
        #     for index_2 in index_lst:
        #         smiles_2_dict = str_to_dict(many_df.at[index_2, pharm_lig_smiles_col], return_str=True)
        #         smiles_2 = smiles_2_dict[many_df.at[index_2, pharm_lig_col]][0]
        #         chem_df.at[index_1, index_2] = get_lig_simi(smiles_1, smiles_2)
        #     compare_bar.progress((i + 1)/len(index_lst))



        with st.expander("Entries Table", expanded=False):

            st.markdown(f"##### Entries Table")

            show_df = rename_st_cols(many_df)

            del show_df[rename_col_dict[pdb_id_col]]
            del show_df[rename_col_dict[modelid_col]]

            show_st_dataframe(show_df)

            entries_file_name = st.text_input(
                label="Entries File Name",
                value=f"{rascore_str}_{entry_table_file}",
            )
            download_st_df(show_df, entries_file_name, "Download Entries Table")            

    write_st_end()