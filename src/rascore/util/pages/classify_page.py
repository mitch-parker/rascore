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

import streamlit as st

from rascore.util.functions.path import get_dir_path

from ..functions.path import (
    delete_path,
    get_file_path,
    get_neighbor_path,
    load_table,
    pages_str,
    data_str,
    rascore_str,
    classify_str,
)
from ..functions.file import result_table_file
from ..functions.gui import (
    get_st_file_path,
    save_st_file,
    show_st_dataframe,
    download_st_df,
    rename_st_cols,
)
from ..pipelines.classify_rascore import classify_rascore


def classify_page():

    st_file_lst = st.file_uploader(
        "Upload RAS Structure(s)", accept_multiple_files=True
    )

    with st.form(key="Classify Form"):
        out_file = st.text_input(label="Classify File Name")
        submit_files = st.form_submit_button(label="Classify Structure(s)")

    if submit_files:
        if len(st_file_lst) > 0:

            try:
                with st.spinner(text="Classifying Structure(s)"):

                    for st_file in st_file_lst:
                        save_st_file(st_file)

                    coord_path_lst = [get_st_file_path(x) for x in st_file_lst]

                    classify_path = get_dir_path(
                        dir_str=f"{rascore_str}_{classify_str}",
                        dir_path=get_neighbor_path(__file__, pages_str, data_str),
                    )

                    classify_rascore(coord_path_lst, out_path=classify_path)

                    st.success("Classified Structure(s)")

                    result_file_path = get_file_path(
                        result_table_file, dir_path=classify_path
                    )

                    df = load_table(result_file_path)

                    df = rename_st_cols(df)

                    show_st_dataframe(df)

                    download_st_df(df, "Download Classification Table", out_file)

                    delete_path(classify_path)
                    for coord_path in coord_path_lst:
                        delete_path(coord_path)

            except:
                st.error("Error Analyzing Inputted Structure(s)")
