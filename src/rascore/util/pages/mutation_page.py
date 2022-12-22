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

from ..scripts.build_dih_table import build_dih_table
from ..scripts.make_facet_plot import make_facet_plot
from ..functions.lig import lig_lst_dict
from ..functions.color import red_hex, gray_hex
from ..functions.table import ( build_col_count_dict, mask_equal, lst_col,
                            extract_int, fix_col, make_dict)
from ..functions.col import (pocket_class_col, pdb_id_col, pdb_code_col,
                            chainid_col, bound_prot_swiss_id_col, pharm_class_col, bound_prot_chainid_col,
                            nuc_class_col, mut_status_col, prot_class_col, interf_class_col,
                            bio_lig_col, ion_lig_col, mem_lig_col, pharm_lig_col, match_class_col, 
                            modelid_col, gene_class_col, prot_class_col,
                            phi_col, psi_col, chi1_col, chi2_col,
                            pharm_lig_col, rename_col_dict)
from ..functions.file import entry_table_file, dih_json_file
from ..functions.path import rascore_str
from ..constants.pharm import pharm_color_dict, sp2_name, other_pharm_name, none_pharm_name
from ..constants.prot import prot_color_dict, none_prot_name, other_prot_name
from ..constants.conf import loop_resid_dict, loop_color_dict, sw1_name, sw2_name, y32_name, y71_name, sw1_color, sw2_color
from ..functions.gui import (load_st_table, show_st_dataframe, write_st_end,
                            show_st_structure, show_st_table, get_html_text, 
                            download_st_df, rename_st_cols, reorder_st_cols)
from ..functions.lst import res_to_lst, lst_nums, str_to_lst

import streamlit as st

reverse_col_dict = make_dict(list(rename_col_dict.values()),list(rename_col_dict.keys()))

sw1_resid_lst = res_to_lst(loop_resid_dict[sw1_name])
sw2_resid_lst = res_to_lst(loop_resid_dict[sw2_name])

def mutation_page():
    
    st.markdown("## Analyze Mutations")

    st.markdown("---")

    df = load_st_table(__file__)

    st.sidebar.markdown("""**Note.** We recommend comparing mutated structures with 
    identical SW1 and SW2 conformations and molecular annotations.
    """)

    left_name = "Left"
    right_name = "Right"

    a_name = "G12D"
    b_name = "G12V"

    left_color = red_hex
    right_color = gray_hex

    left_query_col, right_query_col = st.columns(2)

    mut_status_lst = lst_col(df, mut_status_col, unique=True)

    left_name = left_query_col.selectbox(f"{rename_col_dict[mut_status_col]} ({left_name})", mut_status_lst, index=mut_status_lst.index(a_name))

    mut_status_lst = [x for x in mut_status_lst if x != left_name]

    right_name = right_query_col.selectbox(f"{rename_col_dict[mut_status_col]} ({right_name})", mut_status_lst, index=mut_status_lst.index(b_name))

    query_col_dict = {left_name: left_query_col, right_name: right_query_col}

    query_df_dict = dict()

    for query_name, query_col in query_col_dict.items():

        query_df = df.copy(deep=True)

        query_df = mask_equal(query_df, mut_status_col, query_name)

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

        query_col.success(f"{query_name} (N={len(query_df)})")

    mut_df = pd.concat([query_df_dict[left_name], query_df_dict[right_name]], sort=False)

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
            gene_class = pdb_df[gene_class_col].iloc[0]      

            mut_status = pdb_df[mut_status_col].iloc[0]

            info_col.markdown(f"##### PDB: [{pdb_code.upper()}](https://www.rcsb.org/structure/{pdb_code}) (Chain {chainid}) - {gene_class}({mut_status})")

            col_pair_dict = {nuc_class_col: bio_lig_col, match_class_col: pharm_lig_col, prot_class_col: bound_prot_swiss_id_col}
            for col_1, col_2 in col_pair_dict.items():
                val_1 = pdb_df[col_1].iloc[0]
                val_2 = pdb_df[col_2].iloc[0]
                col_str = f"**{rename_col_dict[col_1]}:** {val_1}"
                if val_1 != "None":
                    if col_2 in [bio_lig_col, pharm_lig_col] and val_2 != "None":
                        col_str += f" ([{val_2}](https://www.rcsb.org/ligand/{val_2}))"
                    elif col_2 in [bound_prot_swiss_id_col] and val_2 != "None":
                        col_str += f" ([{val_2}](https://www.uniprot.org/uniprot/{val_2}))"
                    else:
                        col_str += f" ({val_2})"
                info_col.markdown(col_str)
           
        st.markdown("---")

        st.markdown("##### Viewer Settings") 

        left_set_col, right_set_col = st.columns(2)  

        left_color = left_set_col.color_picker(label=f"{left_name} Color", value=left_color)
        right_color = left_set_col.color_picker(label=f"{right_name} Color", value=right_color)

        stick_resid_lst = left_set_col.multiselect("Displayed Residues (Mutations Always Shown)", lst_nums(1, 189), default=[12, 13, 32, 61])

        label_muts = left_set_col.checkbox("Label Mutations",value=True)
        label_resids = left_set_col.checkbox("Label Residues",value=False)

        style_dict = {"Ribbon": "oval", "Trace": "trace"}

        cartoon_style = style_dict[
            right_set_col.radio("Cartoon Style", ["Ribbon", "Trace"])
        ]

        cartoon_trans = right_set_col.slider(
            "Cartoon Transparency", min_value=0.0, max_value=1.0, value=0.5
        )

        surface_trans = right_set_col.slider(
            "Surface Transparency (Mutations Only)", min_value=0.0, max_value=1.0, value=0.5
        )

        left_view_col, right_view_col = st.columns(2)   

        view_col_dict = {left_name: left_view_col, right_name: right_view_col}  

        for query_name, query_df in query_df_dict.items():
            
            pdb_df = pdb_df_dict[query_name] 

            pdb_code = pdb_df[pdb_code_col].iloc[0]
            chainid = pdb_df[chainid_col].iloc[0]

            mut_status = pdb_df[mut_status_col].iloc[0]
            
            mut_resid_lst = list()
            if mut_status != "WT":
                mut_resid_lst = [extract_int(x) for x in str_to_lst(mut_status)]

            if mut_status == left_name:
                mut_color = left_color
            elif mut_status == right_name:
                mut_color = right_color

            view_col = view_col_dict[query_name]

            show_st_structure(pdb_df,
                    mut_resids=mut_resid_lst,
                    stick_resids=stick_resid_lst,
                    label_resids=label_resids, 
                    label_muts=label_muts,
                    mut_color=mut_color,
                    cartoon_style=cartoon_style,
                    cartoon_trans=cartoon_trans, 
                    mut_trans=surface_trans,
                    zoom_resids=mut_resid_lst + stick_resid_lst,
                    zoom=0.7,
                    width=400,
                    height=400,
                    st_col=view_col)
            
            label_size = "medium"

            for col in [sw1_name, sw2_name, y32_name, y71_name]:

                label_str = get_html_text({f"{rename_col_dict[col]}: ": "#31333F"}, font_weight='bold', font_size=label_size)

                if col in [y32_name, sw1_name]:
                    label_color = loop_color_dict[sw1_name]
                elif col in [y71_name, sw2_name]:
                    label_color = loop_color_dict[sw2_name]

                label_str += get_html_text({pdb_df[col].iloc[0]: label_color}, font_size=label_size)

                view_col.markdown(label_str, unsafe_allow_html=True)

    with st.expander("Many-to-Many Comparison", expanded=False):

        st.markdown("#### Many-to-Many")

        count_dict = build_col_count_dict(mut_df, mut_status_col)

        left_plot_col, middle_plot_col, right_plot_col = st.columns(3)

        left_color = left_plot_col.color_picker(label=f"{left_name} (N={count_dict[left_name]})", value=left_color)
        right_color = left_plot_col.color_picker(label=f"{right_name} (N={count_dict[right_name]})", value=right_color)

        marker_size = left_plot_col.number_input("Marker Size", min_value=1, max_value=50, value=5)

        bb_name = "Backbone Dihedrals (Phi vs. Psi)"
        sc_name = "Side Chain Dihedrals (Chi1 vs. Chi2)"

        plot_type = middle_plot_col.radio("Plot Type", [bb_name, sc_name])

        dih_resids = middle_plot_col.selectbox("Dihedral Residue", lst_nums(1,166), index=31)

        bb_resids = None
        chi1_resids = None
        chi2_resids = None

        if plot_type == bb_name:
            bb_resids = dih_resids
            x_col = phi_col
            y_col = psi_col
        elif plot_type == sc_name:
            chi1_resids = dih_resids
            chi2_resids = dih_resids
            x_col = chi1_col
            y_col = chi2_col

        if middle_plot_col.button("Plot Dihedrals"):

            with st.spinner("Preparing Dihedrals"):

                dih_dict = load_st_table(__file__, file_name=dih_json_file, json_format=True)

                mut_df = build_dih_table(mut_df, dih_dict, bb_resids=bb_resids, 
                                            chi1_resids=chi1_resids,
                                            chi2_resids=chi2_resids)

            fig = make_facet_plot(
                        mut_df,
                        x_col=x_col,
                        y_col=y_col,
                        show_legend=False,
                        hue_col=mut_status_col,
                        hue_order=[left_name, right_name],
                        hue_palette={left_name: left_color, right_name: right_color},
                        plot_width=3,
                        plot_height=3,
                        font_size=12,
                        marker_size=marker_size,
                    )

            right_plot_col.pyplot(fig)

        st.markdown("---")

        left_conf_col, right_conf_col = st.columns(2)

        for table_col in [y32_name, y71_name, sw1_name, sw2_name]:

            if table_col in [y32_name, sw1_name]:
                conf_col = left_conf_col
            elif table_col in [y71_name, sw2_name]:
                conf_col = right_conf_col

            conf_col.markdown(f"##### {rename_col_dict[table_col]}")

            loop_df = (
                pd.pivot_table(
                    data=rename_st_cols(mut_df),
                    index=[rename_col_dict[table_col]],
                    columns=rename_col_dict[mut_status_col],
                    values=rename_col_dict[pdb_id_col],
                    aggfunc="nunique",
                    margins=True,
                )
                .fillna("")
            )

            for col in list(loop_df.columns):
                loop_df[col] = loop_df[col].map(str)
                loop_df = fix_col(loop_df, col)

            loop_df = reorder_st_cols(loop_df, table_col, mut_status_col)

            loop_df = loop_df.reset_index()

            show_st_table(loop_df, st_col=conf_col)

        #     st.markdown("---")

        #     left_site_col, right_chem_col = st.columns(2)

        #     pdb_1 = pdb_name_dict[left_name]
        #     pdb_2 = pdb_name_dict[right_name]

        #     pdb_1_df = pdb_df_dict[left_name]
        #     pdb_2_df = pdb_df_dict[right_name]

        #     lig_1 = pdb_1_df[pharm_lig_col].iloc[0]
        #     lig_2 = pdb_2_df[pharm_lig_col].iloc[0]

        #     smiles_dict_1 = str_to_dict(pdb_1_df[pharm_lig_smiles_col].iloc[0])
        #     smiles_dict_2 = str_to_dict(pdb_2_df[pharm_lig_smiles_col].iloc[0])
            
        #     smiles_1 = smiles_dict_1[lig_1][0]
        #     smiles_2 = smiles_dict_2[lig_2][0]

        #     right_chem_col.markdown("#### Chemistry Comparison")

        #     mcs = get_lig_mcs([smiles_1, smiles_2])

        #     right_chem_col.markdown(f"**2D Fingerprint Similarity:** {round(get_lig_simi(smiles_1, smiles_2),2)}")

        #     right_chem_col.markdown(f"**Maximum Common Substructure (MCS):** {get_lig_smiles(mcs)}")

        #     if not right_chem_col.checkbox("Highlight MCS",value=True):
        #         mcs = None
            
        #     right_chem_col.image(draw_lig_plot([smiles_1, smiles_2], 
        #     lig_labels=[f"{pdb_1} ({lig_1})", f"{pdb_2} ({lig_2})"], font_size=16, mol_pad=0,
        #     plot_height=4, plot_width=3, highlight_querys=mcs, color_palette=[pharm_color_dict[pharm_class]]))

        #     cont_dict_1 = str_to_dict(pdb_1_df[bound_lig_cont_col].iloc[0])
        #     cont_dict_2 = str_to_dict(pdb_2_df[bound_lig_cont_col].iloc[0])

        #     cont_lst_1 = cont_dict_1[lig_1]
        #     cont_lst_2 = cont_dict_2[lig_2]
            
        #     cont_df = pd.DataFrame()

        #     left_site_col.markdown("#### Site Comparison")

        #     left_site_col.markdown(f"**Residue Contact Similarity (Simpson):** {round(calc_simpson(cont_lst_1, cont_lst_2),2)}")

        #     cont_col_1 = f"{pdb_1} ({lig_1})"
        #     cont_col_2 = f"{pdb_2} ({lig_2})"

        #     for resid in sorted(lst_unique(cont_lst_1 + cont_lst_2,return_int=True)):

        #         cont_1 = ''
        #         cont_2 = ''

        #         if str(resid) in cont_lst_1:
        #             cont_1 = "✓"
        #         if str(resid) in cont_lst_2:
        #             cont_2 = "✓"

        #         cont_df.at[resid, cont_col_1] = cont_1
        #         cont_df.at[resid, cont_col_2] = cont_2

        #     cont_df = cont_df.rename_axis('Residue Contact').reset_index()

        #     show_st_table(cont_df, st_col=left_site_col)


        # with st.expander("Many-to-Many Comparison", expanded=False):

        #     st.markdown("#### Many-to-Many Comparison")

        #     count_dict = build_col_count_dict(many_df, id_col)

        #     color_dict = {left_name: pharm_color_dict[pharm_class], 
        #     right_name: get_hex(change_hex_alpha(pharm_color_dict[pharm_class], 0.5)), both_name: gray_hex}

        #     left_plot_col, middle_plot_col, right_plot_col = st.columns(3)

        #     for label_name, label_color in color_dict.items():

        #         count = 0
        #         if label_name in list(count_dict.keys()):
        #             count = count_dict[label_name]

        #         color_dict[label_name] = left_plot_col.color_picker(label=f"{label_name} (N={count})", 
        #                                                         value=label_color)

        #     many_order =  [x for x in list(color_dict.keys()) if x in lst_col(many_df, id_col, unique=True)]

        #     scatter_name = "Scatterplot"
        #     box_name = "Boxplot"

        #     plot_type = middle_plot_col.radio("Plot Type", [scatter_name, box_name])

        #     marker_size = middle_plot_col.number_input("Marker Size", min_value=1, max_value=50, value=5)
            
        #     if plot_type == scatter_name:
        #         plot_reg = middle_plot_col.checkbox("Display Regression",value=True)

        #         fig = make_facet_plot(many_df, 
        #             x_col=pocket_volume_col,
        #             y_col=pocket_score_col,
        #             x_str=rename_col_dict[pocket_volume_col],
        #             y_str=rename_col_dict[pocket_score_col],
        #             x_round=0,
        #             y_round=1,
        #             hue_col=id_col,
        #             hue_palette=color_dict,
        #             hue_order=many_order,
        #             plot_width=3,
        #             plot_height=3,
        #             font_size=12,
        #             x_rotation=45,
        #             marker_size=marker_size,
        #             x_pad=25,
        #             y_pad=10,
        #             x_ha="right",
        #             plot_reg=plot_reg,
        #             log_reg=True,
        #             trun_reg=True,
        #             x_lim=[0, 1700],
        #             x_ticks=[0, 500, 1000, 1500],
        #             y_lim=[0, 1.2],
        #             y_ticks=[0.0, 0.5, 1.0])
            
        #     elif plot_type == box_name:

        #         y_str = middle_plot_col.radio("Y-Axis",[rename_col_dict[pocket_score_col],rename_col_dict[pocket_volume_col]])

        #         y_col = reverse_col_dict[y_str]

        #         if y_col == pocket_score_col:
        #             y_lim = [0, 1.2]
        #             y_ticks = [0.0, 0.5, 1.0]
        #             y_round = 1

        #         elif y_col == pocket_volume_col:
        #             y_lim = [0, 1700]
        #             y_ticks = [0, 500, 1000, 1500]
        #             y_round = 0

        #         stat_pairs = None
        #         if len(many_order) > 1:
        #             if middle_plot_col.checkbox("Display Welch's t-test p-value",value=True):
        #                 stat_pairs = list()
        #                 for id_1 in many_order:
        #                     for id_2 in many_order:
        #                         if id_1 != id_2:
        #                             if (id_2, id_1) not in stat_pairs:
        #                                 stat_pairs.append((id_1, id_2))

        #         fig = make_facet_plot(
        #                 many_df,
        #                 x_col=id_col,
        #                 x_order=many_order,
        #                 x_palette=color_dict,
        #                 x_count=True,
        #                 x_str="",
        #                 y_col=y_col,
        #                 y_str=y_str,
        #                 show_legend=False,
        #                 plot_width=3,
        #                 plot_height=3,
        #                 plot_kind="box",
        #                 font_size=12,
        #                 marker_size=marker_size,
        #                 x_pad=35,
        #                 y_lim=y_lim,
        #                 y_ticks=y_ticks,
        #                 y_round=y_round,
        #                 stat_pairs=stat_pairs,
        #                 stat_test="t-test_welch",
        #                 stat_format="simple"
                
                    
        #             )

        #     right_plot_col.pyplot(fig)

        #     many_df[pocket_score_col] = many_df[pocket_score_col].map(float)
        #     many_df[pocket_volume_col] = many_df[pocket_volume_col].map(float)

        #     left_stat_col, right_stat_col = st.columns(2)

        #     stat_col_dict = {pocket_score_col: left_stat_col, pocket_volume_col: right_stat_col}

        #     for stat_name, stat_col in stat_col_dict.items():

        #         stat_col.markdown(f"##### {rename_col_dict[stat_name]}")

        #         temp_df = pd.pivot_table(many_df, values=stat_name,
        #                                 index=id_col,
        #                                 aggfunc=[np.mean, min, max])

        #         stat_df = pd.DataFrame()
        #         for index in list(temp_df.index.values):
        #             for col, stat in zip(list(temp_df.columns), [mean_name, min_name, max_name]):
        #                 stat_df.at[index, stat] = temp_df.at[index, col]

        #         for col in list(stat_df.columns):
        #             stat_df[col] = stat_df[col].round(1)
        #             stat_df[col] = stat_df[col].map(str)
        #             stat_df = fix_col(stat_df, col)

        #         stat_df = stat_df.loc[many_order, :]

        #         show_st_table(stat_df.rename_axis(rename_col_dict[id_col]).reset_index(), st_col=stat_col)

        #     st.markdown("---")

        #     cont_df = pd.DataFrame()

        #     cont_order = list()

        #     for index in list(many_df.index.values):
                
        #         lig = many_df.at[index, pharm_lig_col]
        #         cont_dict = str_to_dict(many_df.at[index, bound_lig_cont_col], return_int=True)

        #         cont_lst = cont_dict[lig]

        #         label = many_df.at[index, id_col]
        #         label_name = f"{label} (N={count_dict[label]})"

        #         if label_name not in list(cont_df.columns):
        #             cont_df[label_name] = 0

        #         for cont in cont_lst:
        #             if cont in cont_order:
        #                 cont_df.at[cont, label_name] += 1
        #             else:
        #                 cont_df.at[cont, label_name] = 1
        #                 cont_order.append(cont)

        #     cont_df = cont_df.fillna(0)

        #     for col in list(cont_df.columns):
        #         cont_df[col] = (cont_df[col] / count_dict[col.split(" (")[0]]) * 100
        #         cont_df[col] = cont_df[col].round(1)
        #         cont_df[col] = cont_df[col].map(str)
        #         cont_df = fix_col(cont_df, col)
        #         cont_df[col] += "%"
        #         cont_df[col] = cont_df[col].replace("0%","")

        #     cont_df = cont_df.loc[sorted(cont_order), [f"{x} (N={count_dict[x]})" for x in many_order if x in list(count_dict.keys())]]

        #     cont_df = cont_df.rename_axis('Residue Contact').reset_index()

        #     left_table_col, right_table_col = st.columns(2)

        #     right_table_col.markdown("##### Site Comparison")

        #     show_st_table(cont_df, st_col=right_table_col)

        #     for table_col in [y32_name, y71_name, sw1_name, sw2_name]:

        #         left_table_col.markdown(f"##### {rename_col_dict[table_col]}")

        #         loop_df = (
        #             pd.pivot_table(
        #                 data=rename_st_cols(many_df),
        #                 index=[rename_col_dict[table_col]],
        #                 columns=rename_col_dict[id_col],
        #                 values=rename_col_dict[pdb_id_col],
        #                 aggfunc="nunique",
        #                 margins=True,
        #             )
        #             .fillna("")
        #         )

        #         for col in list(loop_df.columns):
        #             loop_df[col] = loop_df[col].map(str)
        #             loop_df = fix_col(loop_df, col)

        #         loop_df = reorder_st_cols(loop_df, table_col, id_col)

        #         loop_df = loop_df.reset_index()

        #         show_st_table(loop_df, st_col=left_table_col)

       
    with st.expander("Entries Table", expanded=False):

        st.markdown(f"##### Entries Table")

        show_df = rename_st_cols(mut_df)

        del show_df[rename_col_dict[pdb_id_col]]
        del show_df[rename_col_dict[modelid_col]]

        show_st_dataframe(show_df)

        entries_file_name = st.text_input(
            label="Entries File Name",
            value=f"{rascore_str}_{entry_table_file}",
        )
        download_st_df(show_df, entries_file_name, "Download Entries Table")            

    write_st_end()