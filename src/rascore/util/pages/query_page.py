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

import pandas as pd
import streamlit as st

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
)
from ..constants.nuc import nuc_class_lst
from ..constants.gene import gene_class_lst
from ..functions.table import lst_col, fix_col, mask_equal
from ..constants.pml import sup_resids, sup_pdb_code, sup_chainid, mono_view
from ..functions.col import (
    rename_col_dict,
    pdb_id_col,
    sw1_col,
    sw2_col,
    mut_status_col,
    prot_class_col,
    pharm_class_col,
    match_class_col,
    pocket_class_col,
    interf_class_col,
    gene_class_col,
    nuc_class_col,
    pdb_code_col,
    core_path_col,
    bio_lig_col,
    ion_lig_col,
    pharm_lig_col,
    chem_lig_col,
    mod_lig_col,
    mem_lig_col,
    pocket_lig_col,
    bound_prot_chainid_col,
    pocket_path_col,
    interf_path_col,
    bound_interf_chainid_col,
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

from ..constants.conf import sw1_name, sw2_name, loop_resid_dict, conf_color_dict
from ..functions.file import (
    entry_table_file,
    sum_table_file,
    pymol_pml_file,
)
from ..scripts.write_pymol_script import pymol_color_dict


def query_page():

    st.markdown("# Query Database")

    st.markdown("---")

    df = load_st_table(__file__)

    st.sidebar.markdown("## Query Selection")

    annot_col_lst = [
        sw1_col,
        sw2_col,
        mut_status_col,
        prot_class_col,
        pharm_class_col,
        match_class_col,
        pocket_class_col,
        interf_class_col,
    ]

    mask_dict = dict()

    mask_df = df.copy(deep=True)

    for col in annot_col_lst:

        mask_dict[col] = st.sidebar.multiselect(
            rename_col_dict[col],
            lst_col(mask_df, col, unique=True),
        )

        mask_df = mask_st_table(mask_df, mask_dict)

    if len(mask_df) == 0:
        st.warning("No Structures Available Based On Inputted Query")
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

        left_table_col, right_table_col = st.columns(2)

        for loop_name in [sw1_name, sw2_name]:

            if loop_name == sw1_name:
                table_col = left_table_col
            elif loop_name == sw2_name:
                table_col = right_table_col

            table_col.markdown(f"#### {loop_name} Conformation(s)")

            loop_df = (
                pd.pivot_table(
                    data=rename_st_cols(gene_nuc_df),
                    index=rename_col_dict[loop_name],
                    columns=rename_col_dict[gene_class_col],
                    values=rename_col_dict[pdb_id_col],
                    aggfunc="nunique",
                    margins=True,
                )
                .reset_index()
                .fillna("")
            )

            for col in list(loop_df.columns):
                loop_df[col] = loop_df[col].map(str)
                loop_df = fix_col(loop_df, col)

            show_st_table(loop_df, st_col=table_col)

        left_table_col.markdown(
            '*Note.* "3P" for GTP or GTP analog-bound, "2P" for GDP-bound, and "0P" for nucleotide-free'
        )

        st.markdown("---")

        st.markdown("#### Download PyMOL Script(s)")

        pymol_lst = [x for x in list(pymol_color_dict.keys()) if x != pocket_lig_col]

        pymol_df = gene_nuc_df.copy(deep=True)

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

        left_pymol_col, right_pymol_col = st.columns(2)

        show_lst = left_pymol_col.multiselect(
            "Show Molecular Contents",
            [rename_col_dict[x] for x in pymol_lst],
        )

        show_color_dict = dict()
        for col in list(pymol_color_dict.keys()):
            rename_col = rename_col_dict[col]
            color = None
            if rename_col in show_lst:
                color = left_pymol_col.text_input(
                    f"{rename_col} Color", value=pymol_color_dict[col]
                )
            show_color_dict[col] = color

        fetch_path = left_pymol_col.text_input(
            label="Fetch Path (e.g., /Users/mitch-parker/rascore)",
        )

        for loop_name, loop_resids in loop_resid_dict.items():

            if loop_name == sw1_name:
                stick_resid = [32]
            elif loop_name == sw2_name:
                stick_resid = [71]

            pymol_file_path = get_file_path(
                f"{loop_name}_{pymol_pml_file}",
                dir_path=get_neighbor_path(__file__, pages_str, data_str),
            )

            pymol_file_name = right_pymol_col.text_input(
                label=f"{loop_name} Script File Name",
                value=f"{loop_name}_{pymol_pml_file}",
            )
            if right_pymol_col.button(f"Create {loop_name} Script File"):
                with st.spinner(text=f"Creating {loop_name} PyMOL File"):
                    write_pymol_script(
                        pymol_df,
                        pymol_file_path,
                        group_col=loop_name,
                        stick_resids=stick_resid,
                        loop_resids=loop_resids,
                        style_ribbon=True,
                        thick_bb=False,
                        color_palette=conf_color_dict[loop_name],
                        show_bio=show_color_dict[bio_lig_col],
                        show_ion=show_color_dict[ion_lig_col],
                        show_pharm=show_color_dict[pharm_lig_col],
                        show_chem=show_color_dict[chem_lig_col],
                        show_mod=show_color_dict[mod_lig_col],
                        show_mem=show_color_dict[mem_lig_col],
                        show_pocket=show_color_dict[pocket_lig_col],
                        show_prot=show_color_dict[bound_prot_chainid_col],
                        sup_resids=sup_resids,
                        show_resids=sup_resids,
                        sup_coord_path=sup_coord_path,
                        sup_chainid=sup_chainid,
                        set_view=mono_view,
                        coord_path_col=coord_path_col,
                        fetch_path=fetch_path,
                    )

                download_st_file(
                    pymol_file_path,
                    pymol_file_name,
                    f"Download {loop_name} Script File",
                    st_col=right_pymol_col,
                )

                delete_path(pymol_file_path)

        st.markdown("---")

        st.markdown("#### Summary Table")

        left_sum_col, right_sum_col = st.columns(2)

        sum_col_lst = [gene_class_col, nuc_class_col] + annot_col_lst

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
                    .reset_index()
                    .fillna("")
                )

                for col in list(sum_df.columns):
                    sum_df[col] = sum_df[col].map(str)
                    sum_df = fix_col(sum_df, col)

                show_st_table(sum_df)

                sum_file_name = st.text_input(
                    label="Summary File Name", value=f"{rascore_str}_{sum_table_file}"
                )
                download_st_df(sum_df, sum_file_name, "Download Summary Table")

                st.markdown("---")

                st.markdown("#### Entries Table")

                show_st_dataframe(gene_nuc_df)

                entries_file_name = st.text_input(
                    label="Entries File Name", value=f"{rascore_str}_{entry_table_file}"
                )
                download_st_df(gene_nuc_df, entries_file_name, "Download Entries Table")

        write_st_end()
