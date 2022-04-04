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
from random import randint

from ..functions.col import id_col, core_path_col, pharm_lig_site_col, pocket_class_col
from ..functions.path import (
    delete_path,
    get_file_path,
    get_file_name,
    get_neighbor_path,
    load_table,
    get_dir_path,
    pages_str,
    data_str,
    rascore_str,
    classify_str,
)
from ..functions.lst import res_to_lst
from ..functions.file import result_table_file
from ..functions.gui import (
    save_st_file,
    show_st_dataframe,
    download_st_df,
    rename_st_cols,
    write_st_end,
)
from ..pipelines.classify_rascore import classify_rascore


def classify_page():

    st.markdown("## Classify Structures")

    st.markdown("---")

    st_file_lst = st.file_uploader(
        "Upload Coordinate Files", accept_multiple_files=True
    )

    table_st_file = None
    with st.expander("Optional Input", expanded=False):

        left_col, right_col = st.columns(2)

        sw1_resids = left_col.text_input('SW1 Residues',value="25-40")
        sw2_resids = right_col.text_input('SW2 Residues',value="56-76")

        sw1_len = len(res_to_lst('25-40'))
        sw2_len = len(res_to_lst('57-76'))

        if len(res_to_lst(sw1_resids)) != sw1_len:
            left_col.warning(f'SW1 Must Be {sw1_len} Residue Long')

        if len(res_to_lst(sw2_resids)) != len(res_to_lst('56-76')):
            right_col.warning(f'SW2 Must Be {sw2_len} Residue Long')

        y32_resid = int(left_col.text_input('Y32 Residue',value=32))
        y71_resid = int(right_col.text_input('Y71 Residue',value=71))

        g12_resid = int(left_col.text_input('G12 Residue',value=12))
        v9_resid = int(right_col.text_input('V9 Residue',value=9))

        st.markdown(
            """
        Tab-separated table with columns:
        - **core_path:** coordinate path:
        - **modelid:** model number (*optional*)
        - **chainid:** chain identifier
        - **nucleotide_class:** nucleotide state (*optional*)
        """
        )
        table_st_file = st.file_uploader(
            "Upload Table File", accept_multiple_files=False
        )

    with st.form(key="Classify Form"):
        out_file = st.text_input(
            label="Classify File Name", value="rascore_classify.tsv"
        )
        submit_files = st.form_submit_button(label="Classify Structures")

    if submit_files:
        if len(st_file_lst) > 0:
            with st.spinner(text="Classifying Structures"):

                try:
                    coord_path_lst = list()
                    id_dict = dict()
                    path_dict = dict()
                    for st_file in st_file_lst:
                        coord_path = save_st_file(st_file)
                        coord_path_lst.append(coord_path)
                        id_dict[get_file_name(coord_path)] = st_file.name
                        path_dict[st_file.name] = coord_path

                    classify_path = get_dir_path(
                        dir_str=f"{rascore_str}_{classify_str}_{randint(0,3261994)}",
                        dir_path=get_neighbor_path(__file__, pages_str, data_str),
                    )

                    if table_st_file is None:
                        classify_input = coord_path_lst
                    else:
                        table_path = save_st_file(table_st_file)
                        classify_input = load_table(table_path)
                        classify_input[core_path_col] = classify_input[
                            core_path_col
                        ].map(path_dict)
                        delete_path(table_path)

                    classify_rascore(classify_input, y32_resid=y32_resid, y71_resid=y71_resid, 
                    g12_resid=g12_resid,v9_resid=v9_resid,
                    sw1_resids=sw1_resids, sw2_resids=sw2_resids, 
                    out_path=classify_path)

                    st.success("Classified Structures")

                    result_file_path = get_file_path(
                        result_table_file, dir_path=classify_path
                    )

                    df = load_table(result_file_path)

                    df[id_col] = df[id_col].map(id_dict)

                    df = df.rename(columns={pharm_lig_site_col: pocket_class_col})

                    df = rename_st_cols(df)

                    show_st_dataframe(df)

                    download_st_df(df, out_file, "Download Classification Table")
                except:
                    st.error("Error Analyzing Uploaded Structures")

                delete_path(classify_path)
                for coord_path in coord_path_lst:
                    delete_path(coord_path)

    write_st_end()
