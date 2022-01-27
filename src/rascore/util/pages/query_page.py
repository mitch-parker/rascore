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
    download_st_text,
    show_st_dataframe,
    write_st_end,
)
from ..constants.gene import gene_class_lst
from ..functions.table import lst_col, fix_col
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
)
from ..functions.path import get_file_path, get_neighbor_path, pages_str, data_str
from ..constants.conf import sw1_name, sw2_name, sw1_resids, sw2_resids, conf_color_dict
from ..functions.file import pymol_pml_file


def query_page():

    st.markdown("# Query Database")

    df = load_st_table(__file__)
    df = rename_st_cols(df)

    st.sidebar.markdown("## Query Selection")

    annot_lst = [
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

    for annot in annot_lst:

        mask_dict[rename_col_dict[annot]] = st.sidebar.multiselect(
            rename_col_dict[annot],
            lst_col(mask_df, rename_col_dict[annot], unique=True),
        )

        mask_df = mask_st_table(mask_df, mask_dict)

    if len(mask_df) == 0:
        st.warning("No Structures Available Based On Inputted Query")
    else:

        left_query_col, right_query_col = st.columns(2)

        gene_lst = lst_col(mask_df, rename_col_dict[gene_class_col])

        radio_lst = [f"All (N={len(gene_lst)})"]
        for gene in gene_class_lst:
            total_gene = len([x for x in gene_lst if x == gene])
            if total_gene > 0:
                radio_lst.append(f"{gene} (N={total_gene})")

        gene_name = left_query_col.radio("Select RAS Structures", radio_lst)

        gene_df = mask_st_table(
            mask_df, {rename_col_dict[gene_class_col]: gene_name.split(" (")[0]}
        )

        right_query_col.markdown("#### Download PyMOL Scripts")

        with right_query_col.form(key="SW1 Form"):
            sw1_file_name = st.text_input(label="SW1 PyMOL File Name")
            sw1_submit = st.form_submit_button(label="Create SW1 PyMOL File")

        if sw1_submit:
            sw1_pymol_path = get_file_path(
                f"{sw1_name}_{pymol_pml_file}",
                dir_path=get_neighbor_path(__file__, pages_str, data_str),
            )
            with st.spinner(text="Creating SW1 PyMOL File"):
                write_pymol_script(gene_df, sw1_pymol_path)

        with right_query_col.form(key="SW2 Form"):
            sw2_file_name = st.text_input(label="SW2 PyMOL File Name")
            sw2_submit = st.form_submit_button(label="Create SW2 PyMOL File")

        if sw2_submit:
            sw2_pymol_path = get_file_path(
                f"{sw2_name}_{pymol_pml_file}",
                dir_path=get_neighbor_path(__file__, pages_str, data_str),
            )
            with st.spinner(text="Creating SW2 PyMOL File"):
                write_pymol_script(gene_df, sw2_pymol_path)

        st.markdown("#### Pivot Table")

        left_pivot_col, right_pivot_col = st.columns(2)

        row_lst = left_pivot_col.multiselect(
            "Rows", [rename_col_dict[x] for x in annot_lst]
        )

        col_lst = right_pivot_col.multiselect(
            "Columns",
            [rename_col_dict[x] for x in annot_lst],
        )

        if len(row_lst) == 0 and len(col_lst) == 0:
            col_lst = list([rename_col_dict[sw1_col]])
            row_lst = list([rename_col_dict[sw2_col]])

        if len([x for x in row_lst if x in col_lst]) > 0:
            st.warning("Cannot Have Same Selection in Rows and Columns")
        else:
            if len(row_lst) > 0 and len(col_lst) > 0:
                pivot_df = (
                    pd.pivot_table(
                        data=gene_df,
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

                pivot_file_name = st.text_input(label="Pivot File Name")
                download_st_df(pivot_df, pivot_file_name, "Download Pivot Table")

                st.markdown("---")

                st.markdown("#### Entries Table")

                show_st_dataframe(mask_df)

                entries_file_name = st.text_input(label="Entries File Name")
                download_st_df(gene_df, entries_file_name, "Download Entries Table")

        write_st_end()
