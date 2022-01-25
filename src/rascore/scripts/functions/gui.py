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
import py3Dmol
from stmol import showmol

from .file import entry_table_file
from .table import mask_equal
from .path import (
    get_dir_name,
    get_dir_path,
    get_file_path,
    load_table,
    pages_str,
    data_str,
)
from .col import *
from .lst import *


def load_st_table(file_path):

    dir_path = get_dir_name(file_path)
    dir_path = dir_path.split(pages_str)[0]
    dir_path += data_str

    return load_table(
        get_file_path(
            entry_table_file,
            dir_path=get_dir_path(dir_path=dir_path),
        )
    )


@st.cache
def mask_st_table(df, col_dict):

    mask_df = df.copy()

    for col in col_dict.keys():
        mask_df = mask_equal(mask_df, col, col_dict[col])

    return mask_df


def rename_st_cols(df, col_lst=None):

    if col_lst is None:
        col_lst = list(rename_col_dict.keys())

    return df.loc[:, col_lst].rename(columns=rename_col_dict)


def show_st_dataframe(df, st_col=None):

    hide_dataframe_row_index = """
            <style>
            .row_heading.level0 {display:none}
            .blank {display:none}
            </style>
            """

    st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)

    if st_col is None:
        st.dataframe(df)
    else:
        st_col.dataframe(df)


def show_st_table(df, st_col=None):

    hide_table_row_index = """
            <style>
            tbody th {display:none}
            .blank {display:none}
            </style>
            """
    st.markdown(hide_table_row_index, unsafe_allow_html=True)

    if st_col is None:
        st.table(df)
    else:
        st_col.table(df)


def show_st_structure(
    pdb_code,
    style_lst=None,
    label_lst=None,
    zoom_dict=None,
    surface_lst=None,
    cartoon_style="trace",
    cartoon_radius=0.2,
    cartoon_color="lightgray",
    zoom=1,
    width=700,
    height=500,
):
    view = py3Dmol.view(query=f"pdb:{pdb_code}", width=width, height=height)

    view.setStyle(
        {
            "cartoon": {
                "style": cartoon_style,
                "color": cartoon_color,
                "thickness": cartoon_radius,
            }
        }
    )

    view.setViewStyle({"style": "outline", "color": "black", "width": 0.05})

    if surface_lst is not None:
        for surface in surface_lst:
            view.addSurface(py3Dmol.VDW, surface[0], surface[1])

    if style_lst is not None:
        for style in style_lst:
            view.addStyle(
                style[0],
                style[1],
            )

    if label_lst is not None:
        for label in label_lst:
            view.addLabel(label[0], label[1], label[2])

    if zoom_dict is None:
        view.zoomTo()
    else:
        view.zoomTo(zoom_dict)
    view.zoom(zoom)
    showmol(view, height=height, width=width)
