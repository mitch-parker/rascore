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
from random import randint

from ..scripts.write_pymol_script import write_pymol_script

from ..functions.gui import (
    load_st_table,
    rename_st_cols,
    mask_st_table,
    show_st_table,
    download_st_df,
    download_st_file,
    show_st_dataframe,
    write_st_end,
    reorder_st_cols
)
from ..functions.lst import lst_nums
from ..constants.nuc import nuc_class_lst
from ..constants.gene import gene_class_lst
from ..functions.table import lst_col, fix_col, mask_equal, make_dict
from ..functions.color import get_lst_colors, get_hex
from ..constants.pml import sup_resids, show_resids, sup_pdb_code, sup_chainid, mono_view
from ..functions.col import (
    rename_col_dict,
    pdb_id_col,
    sw1_col,
    sw2_col,
    y32_col,
    y71_col,
    mut_status_col,
    prot_class_col,
    match_class_col,
    pocket_class_col,
    interf_class_col,
    gene_class_col,
    nuc_class_col,
    pdb_code_col,
    modelid_col,
    core_path_col,
    bio_lig_col,
    ion_lig_col,
    pharm_lig_col,
    chem_lig_col,
    mod_lig_col,
    mem_lig_col,
    pocket_lig_col,
    bound_prot_chainid_col,
)
from ..functions.path import (
    get_file_path,
    get_neighbor_path,
    delete_path,
    path_exists,
    pages_str,
    data_str,
    rascore_str,
)

from ..constants.conf import y32_name, y71_name, sw1_name, sw2_name, loop_resid_dict, conf_color_dict, sw1_color,sw2_color
from ..functions.file import (
    entry_table_file,
    sum_table_file,
    pymol_pml_file,
)
from ..scripts.write_pymol_script import pymol_color_dict

reverse_col_dict = make_dict(list(rename_col_dict.values()),list(rename_col_dict.keys()))


def query_page():

    st.markdown("## Query Database")

    st.markdown("---")

    df = load_st_table(__file__)

    conf_col_lst = [sw1_col, sw2_col, y32_col, y71_col]

    annot_col_lst = [
        mut_status_col,
        prot_class_col,
        pocket_class_col,
        match_class_col,
        interf_class_col,
    ]

    mask_dict = dict()

    mask_df = df.copy(deep=True)

    st.sidebar.markdown(
        "**Note.** Selections dynamically update from top to bottom. Multiple selections are possible."
    )

    st.sidebar.markdown("## Conformation Selection")

    for col in conf_col_lst:

        mask_dict[col] = st.sidebar.multiselect(
            rename_col_dict[col],
            lst_col(mask_df, col, unique=True),
        )

        mask_df = mask_st_table(mask_df, mask_dict)

    st.sidebar.markdown("## Annotation Selection")

    for col in annot_col_lst:

        mask_dict[col] = st.sidebar.multiselect(
            rename_col_dict[col],
            lst_col(mask_df, col, unique=True),
        )

        mask_df = mask_st_table(mask_df, mask_dict)

    if len(mask_df) == 0:
        st.warning("No Structures Available Based On Selected Query")
    else:

        left_query_col, right_query_col = st.columns(2)

        gene_lst = lst_col(mask_df, gene_class_col)

        gene_radio_lst = [f"All (N={len(gene_lst)})"]
        for gene in gene_class_lst:
            total_gene = len([x for x in gene_lst if x == gene])
            if total_gene > 0:
                gene_radio_lst.append(f"{gene} (N={total_gene})")

        gene_class = left_query_col.radio(
            f"{rename_col_dict[gene_class_col]}", gene_radio_lst
        )

        gene_df = mask_st_table(mask_df, {gene_class_col: gene_class.split(" (")[0]})

        nuc_lst = lst_col(gene_df, nuc_class_col)

        nuc_radio_lst = [f"All (N={len(nuc_lst)})"]
        for nuc in nuc_class_lst:
            total_nuc = len([x for x in nuc_lst if x == nuc])
            if total_nuc > 0:
                nuc_radio_lst.append(f"{nuc} (N={total_nuc})")

        nuc_class = right_query_col.radio(
            f"{rename_col_dict[nuc_class_col]}", nuc_radio_lst
        )

        gene_nuc_df = mask_st_table(gene_df, {nuc_class_col: nuc_class.split(" (")[0]})

        left_conf_col, right_conf_col = st.columns(2)

        for table_col in [y32_name, y71_name, sw1_name, sw2_name]:

            if table_col in [y32_name, sw1_name]:
                conf_col = left_conf_col
            elif table_col in [y71_name, sw2_name]:
                conf_col = right_conf_col

            conf_col.markdown(f"##### {rename_col_dict[table_col]}")

            loop_df = (
                pd.pivot_table(
                    data=rename_st_cols(gene_nuc_df),
                    index=[rename_col_dict[table_col]],
                    columns=rename_col_dict[gene_class_col],
                    values=rename_col_dict[pdb_id_col],
                    aggfunc="nunique",
                    margins=True,
                )
                .fillna("")
            )

            for col in list(loop_df.columns):
                loop_df[col] = loop_df[col].map(str)
                loop_df = fix_col(loop_df, col)

            loop_df = reorder_st_cols(loop_df, table_col, gene_class_col)

            loop_df = loop_df.reset_index()

            show_st_table(loop_df, st_col=conf_col)

        with st.expander("PyMOL Script", expanded=False):

            st.markdown("#### PyMOL Script")

            pymol_lst = [x for x in list(pymol_color_dict.keys()) if x != pocket_lig_col]

            coord_path_col = pdb_code_col
            sup_coord_path = sup_pdb_code
            if len(
                [x for x in lst_col(gene_nuc_df, core_path_col) if path_exists(x)]
            ) == len(gene_nuc_df):
                coord_path_col = core_path_col
                sup_coord_path = lst_col(
                    mask_equal(df, pdb_id_col, f"{sup_pdb_code}{sup_chainid}"),
                    core_path_col,
                )[0]

            group_col = st.selectbox("Group By Selection", [rename_col_dict[x] for x in conf_col_lst + annot_col_lst])

            show_lst = st.multiselect(
                "Show Molecular Contents",
                [rename_col_dict[x] for x in pymol_lst],default=[rename_col_dict[bio_lig_col],rename_col_dict[pharm_lig_col]]
            )

            left_pymol_col, right_pymol_col = st.columns(2)

            show_color_dict = dict()
            for col in list(pymol_color_dict.keys()):
                rename_col = rename_col_dict[col]
                color = None
                if rename_col in show_lst:
                    color = left_pymol_col.color_picker(f"{rename_col} Color", get_hex(pymol_color_dict[col]))

                show_color_dict[col] = color

            style_dict = {"Ribbon": False, "Trace": True}

            style_ribbon = style_dict[
                right_pymol_col.radio("Cartoon Style", ["Ribbon", "Trace"])
            ]


            stick_resids = right_pymol_col.multiselect(
                'Stick Residues',
                lst_nums(1,189),
                default=[32,71]
            )

            loop_resids = list()
            for loop_name in list(loop_resid_dict.keys()):
                if right_pymol_col.checkbox(f"Highlight {loop_name}",value=True):
                    loop_resids.append(loop_resid_dict[loop_name])
            if len(loop_resids) == 0:
                loop_resids = [show_resids]
                

            color_group = left_pymol_col.checkbox("Color Loops By Group")
            left_pymol_col.write("**Default:** Colors SW1 pink and SW2 purple")

            if color_group:
                if group_col == rename_col_dict[sw1_name]: 
                    color_palette = conf_color_dict[sw1_name]
                elif group_col == rename_col_dict[sw2_name]: 
                    color_palette = conf_color_dict[sw2_name]
                else:
                    group_color_palette = get_lst_colors(lst_col(rename_st_cols(gene_nuc_df),group_col,unique=True),return_dict=True)
                    if len(group_color_palette.items()) > 10:
                        color_palette = group_color_palette.copy()
                    else:
                        left_color_col, right_color_col = st.columns(2)
                        color_palette = dict()
                        i = 0
                        for group, color in group_color_palette.items():
                            if i == 0:
                                color_palette[group] = left_color_col.color_picker(f"{group} Color", get_hex(color))
                                i += 1
                            elif i == 1:
                                color_palette[group] = right_color_col.color_picker(f"{group} Color", get_hex(color))
                                i = 0
            else:
                color_palette = [sw1_color,sw2_color]

            fetch_path = st.text_input(
                label="Fetch Path (e.g., /Users/mitch-parker/rascore)",
            )

            pymol_file_path = get_file_path(
                f"{pymol_pml_file}_{randint(0,3261994)}",
                dir_path=get_neighbor_path(__file__, pages_str, data_str),
            )

            pymol_file_name = st.text_input(
                label="PyMOL Script Name",
                value=f"{rascore_str}_{pymol_pml_file}",
            )
            if st.button("Create PyMOL Script"):
                with st.spinner(text="Creating PyMOL Script"):
                    write_pymol_script(
                        gene_nuc_df,
                        pymol_file_path,
                        group_col=[x for x in conf_col_lst + annot_col_lst if rename_col_dict[x] == group_col][0],
                        stick_resids=stick_resids,
                        loop_resids=loop_resids,
                        style_ribbon=style_ribbon,
                        thick_bb=False,
                        color_group=color_group,
                        color_palette=color_palette,
                        show_bio=show_color_dict[bio_lig_col],
                        show_ion=show_color_dict[ion_lig_col],
                        show_pharm=show_color_dict[pharm_lig_col],
                        show_chem=show_color_dict[chem_lig_col],
                        show_mod=show_color_dict[mod_lig_col],
                        show_mem=show_color_dict[mem_lig_col],
                        show_pocket=show_color_dict[pocket_lig_col],
                        show_prot=show_color_dict[bound_prot_chainid_col],
                        sup_resids=sup_resids,
                        show_resids=show_resids,
                        sup_coord_path=sup_coord_path,
                        sup_chainid=sup_chainid,
                        set_view=mono_view,
                        coord_path_col=coord_path_col,
                        fetch_path=fetch_path,
                    )

                download_st_file(
                    pymol_file_path,
                    pymol_file_name,
                    f"Download PyMOL Script",
                )

                delete_path(pymol_file_path)

        with st.expander("Summary Table", expanded=False):

            st.markdown("#### Summary Table")

            left_sum_col, right_sum_col = st.columns(2)

            sum_col_lst = [gene_class_col, nuc_class_col] + conf_col_lst + annot_col_lst

            row_lst = left_sum_col.multiselect(
                "Rows",
                [rename_col_dict[x] for x in sum_col_lst],
                rename_col_dict[sw2_col],
            )

            col_lst = right_sum_col.multiselect(
                "Columns",
                [rename_col_dict[x] for x in sum_col_lst],
                default=rename_col_dict[sw1_col],
            )

            gene_nuc_df = rename_st_cols(gene_nuc_df)

            if len([x for x in row_lst if x in col_lst]) > 0:
                st.warning("Cannot Have Same Selection in Rows and Columns")
            else:
                if len(row_lst) > 0 and len(col_lst) > 0:
                    sum_df = (
                        pd.pivot_table(
                            data=gene_nuc_df,
                            index=row_lst,
                            columns=col_lst,
                            values=rename_col_dict[pdb_id_col],
                            aggfunc="nunique",
                            margins=True,
                        )
                        .fillna("")
                    )

                    for col in list(sum_df.columns):
                        sum_df[col] = sum_df[col].map(str)
                        sum_df = fix_col(sum_df, col)

                
                    if len(row_lst) == 1 and len(col_lst) == 1:
                        sum_df = reorder_st_cols(sum_df, reverse_col_dict[row_lst[0]], reverse_col_dict[col_lst[0]])

                    sum_df = sum_df.reset_index()

                    show_st_table(sum_df)

                    sum_file_name = st.text_input(
                        label="Summary File Name",
                        value=f"{rascore_str}_{sum_table_file}",
                    )

                    download_st_df(sum_df, sum_file_name, "Download Summary Table")

        with st.expander("Entries Table", expanded=False):

            st.markdown("#### Entries Table")

            del gene_nuc_df[rename_col_dict[pdb_id_col]]
            del gene_nuc_df[rename_col_dict[modelid_col]]

            show_st_dataframe(gene_nuc_df)

            entries_file_name = st.text_input(
                label="Entries File Name",
                value=f"{rascore_str}_{entry_table_file}",
            )
            download_st_df(gene_nuc_df, entries_file_name, "Download Entries Table")

        write_st_end()
