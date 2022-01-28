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
from ..functions.table import lst_col, fix_col, make_dict
from ..constants.pml import sup_resids
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
)
from ..functions.path import (
    get_file_path,
    get_neighbor_path,
    delete_path,
    path_exists,
    pages_str,
    data_str,
)
from ..functions.lig import lig_col_lst
from ..constants.conf import sw1_name, sw2_name, loop_resid_dict, conf_color_dict
from ..functions.file import pymol_pml_file


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

        nuc_class = left_query_col.radio(
            f"{rename_col_dict[nuc_class_col]}", nuc_radio_lst
        )

        left_query_col.markdown(
            '*Note.* "3P" for GTP or GTP analog-bound, "2P" for GDP-bound, and "0P" for nucleotide-free'
        )

        gene_nuc_df = mask_st_table(
            gene_df, {rename_col_dict[nuc_class_col]: nuc_class.split(" (")[0]}
        )

        right_query_col.markdown("#### Download PyMOL Script(s)")

        show_lst = right_query_col.multiselect(
            "Show Molecular Contents",
            [rename_col_dict[x] for x in lig_col_lst],
        )

        for lig_col in lig_col_lst:
            rename_col = rename_col_dict[lig_col]
            if rename_col in show_lst:
                show_bio = right_query_col.text_input(
                    f"{rename_col} Color", value="white"
                )

        right_query_col.markdown("---")

        coord_path_col = pdb_code_col
        if len(
            [x for x in lst_col(gene_nuc_df, core_path_col) if path_exists(x)]
        ) == len(gene_nuc_df):
            coord_path_col = core_path_col

        for loop_name, loop_resids in loop_resid_dict.items():

            if loop_name == sw1_name:
                stick_resid = [32]
            elif loop_name == sw2_name:
                stick_resid = [71]

            pymol_file_path = get_file_path(
                f"{loop_name}_{pymol_pml_file}",
                dir_path=get_neighbor_path(__file__, pages_str, data_str),
            )

            pymol_file_name = right_query_col.text_input(
                label=f"{loop_name} Script File Name",
                value=f"{loop_name}_{pymol_pml_file}",
            )
            if right_query_col.button(f"Create {loop_name} Script File"):
                with st.spinner(text=f"Creating {loop_name} PyMOL File"):
                    write_pymol_script(
                        gene_nuc_df,
                        pymol_file_path,
                        group_col=loop_name,
                        stick_resids=stick_resid,
                        loop_resids=loop_resids,
                        style_ribbon=True,
                        thick_bb=False,
                        color_palette=conf_color_dict[loop_name],
                        show_bio=True,
                        sup_resids=sup_resids,
                        show_resids=sup_resids,
                        coord_path_col=coord_path_col,
                    )

                download_st_file(
                    pymol_file_path,
                    pymol_file_name,
                    f"Download {loop_name} Script File",
                    st_col=right_query_col,
                )

                delete_path(pymol_file_path)

        right_query_col.markdown(
            '*Hint.* When using [PyMOL](https://pymol.org/2/), run the command "set fetch_path [path to coordinate files]" to control the output directory path.'
        )

        st.markdown("---")

        st.markdown("#### Pivot Table")

        left_pivot_col, right_pivot_col = st.columns(2)

        pivot_col_lst = [gene_class_col, nuc_class_col] + annot_col_lst

        row_lst = left_pivot_col.multiselect(
            "Rows",
            [rename_col_dict[x] for x in pivot_col_lst],
            rename_col_dict[sw2_col],
        )

        col_lst = right_pivot_col.multiselect(
            "Columns",
            [rename_col_dict[x] for x in pivot_col_lst],
            default=rename_col_dict[sw1_col],
        )

        gene_nuc_df = rename_st_cols(gene_nuc_df)

        if len([x for x in row_lst if x in col_lst]) > 0:
            st.warning("Cannot Have Same Selection in Rows and Columns")
        else:
            if len(row_lst) > 0 and len(col_lst) > 0:
                pivot_df = (
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

                for col in list(pivot_df.columns):
                    pivot_df[col] = pivot_df[col].map(str)
                    pivot_df = fix_col(pivot_df, col)

                show_st_table(pivot_df)

                pivot_file_name = st.text_input(
                    label="Pivot File Name", value="rascore_pivot.tsv"
                )
                download_st_df(pivot_df, pivot_file_name, "Download Pivot Table")

                st.markdown("---")

                st.markdown("#### Entries Table")

                show_st_dataframe(gene_nuc_df)

                entries_file_name = st.text_input(
                    label="Entries File Name", value="rascore_entries.tsv"
                )
                download_st_df(gene_nuc_df, entries_file_name, "Download Entries Table")

        write_st_end()
