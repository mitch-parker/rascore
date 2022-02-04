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

from ..constants.conf import (
    loop_resid_dict,
    sw1_name,
    sw2_name,
    noise_name,
    conf_color_dict,
)
from ..constants.nuc import nuc_color_dict, nuc_class_lst, gtp_name

from ..functions.gui import (
    get_neighbor_path,
    load_st_table,
    get_html_text,
    write_st_end,
    show_st_structure,
    rename_st_cols,
    show_st_table,
)
from ..functions.path import (
    get_file_path,
    load_table,
    pages_str,
    data_str,
    rascore_str,
    cluster_str,
)
from ..functions.file import result_table_file
from ..functions.col import (
    rename_col_dict,
    nuc_class_col,
    gene_class_col,
    match_class_col,
    prot_class_col,
    pdb_id_col,
    complete_col,
    bio_lig_col,
)
from ..functions.table import (
    mask_equal,
    mask_unequal,
    get_col_most_common,
    lst_col,
    fix_col,
)


def conformation_page():

    df = load_st_table(__file__)

    st.markdown("# Explore Conformations")
    st.markdown("---")

    st.sidebar.markdown(
        """
    **Note.** Only structures with completely modeled loops from original clustering are displayed.
    Tables below structures include updated counts from *rascore* database.
    """
    )

    data_path = get_neighbor_path(__file__, pages_str, data_str)

    sw1_df = load_table(
        get_file_path(
            f"{sw1_name}_{result_table_file}",
            dir_path=f"{data_path}/{rascore_str}_{cluster_str}/{sw1_name}",
        )
    )
    sw2_df = load_table(
        get_file_path(
            f"{sw2_name}_{result_table_file}",
            dir_path=f"{data_path}/{rascore_str}_{cluster_str}/{sw2_name}",
        )
    )
    sw1_df = mask_equal(sw1_df, complete_col, str(True))
    sw2_df = mask_equal(sw2_df, complete_col, str(True))

    table_col_lst = [rename_col_dict[prot_class_col], rename_col_dict[match_class_col]]

    for nuc_class in nuc_class_lst:
        with st.expander(f"{nuc_class} Conformations", expanded=True):
            st.markdown(
                get_html_text(
                    {
                        nuc_class: nuc_color_dict[nuc_class],
                        " Conformations": "#31333F",
                    },
                    font_size="x-large",
                    font_weight="bold",
                ),
                unsafe_allow_html=True,
            )

            sw1_col, sw2_col = st.columns(2)

            nuc_df = mask_equal(df, nuc_class_col, nuc_class)

            for loop_name, loop_resids in loop_resid_dict.items():

                if loop_name == sw1_name:
                    cluster_df = sw1_df.copy(deep=True)
                    stick_resid = 32
                    loop_col = sw1_col

                elif loop_name == sw2_name:
                    cluster_df = sw2_df.copy(deep=True)
                    stick_resid = 71
                    loop_col = sw2_col

                loop_df = mask_unequal(nuc_df, loop_name, noise_name)

                conf_lst = get_col_most_common(loop_df, loop_name)

                if len(conf_lst) > 1:
                    conf_str = f"{loop_name} Conformations ({len(conf_lst)} in Total)"
                else:
                    conf_str = f"{loop_name} Conformation (Only 1)"

                loop_col.markdown(f"#### {conf_str}")

                loop_conf = loop_col.selectbox("Cluster Name", conf_lst)

                pdb_id = loop_col.selectbox(
                    "PDB ID",
                    [
                        x.upper()
                        for x in lst_col(loop_df, pdb_id_col)
                        if x in lst_col(cluster_df, pdb_id_col)
                    ],
                )

                pdb_code = pdb_id[:4].lower()
                chainid = pdb_id[4:5]

                style_lst = list()

                loop_color = conf_color_dict[loop_name][loop_conf]

                cartoon_style = "oval"

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

                bio_lig = lst_col(
                    mask_equal(df, pdb_id_col, f"{pdb_code}{chainid}"), bio_lig_col
                )[0]
                style_lst.append(
                    [
                        {
                            "resn": [bio_lig],
                        },
                        {
                            "stick": {
                                "colorscheme": "whiteCarbon",
                                "radius": 0.2,
                            }
                        },
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
                            }
                        },
                    ]
                )

                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "resi": [stick_resid],
                            "elem": "C",
                        },
                        {"stick": {"color": loop_color, "radius": 0.2}},
                    ]
                )

                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "resi": [stick_resid],
                            "elem": ["O", "N", "H"],
                        },
                        {"stick": {"colorscheme": "Carbon", "radius": 0.2}},
                    ]
                )

                with loop_col:
                    show_st_structure(
                        pdb_code,
                        zoom_dict={"chain": chainid, "resi": [loop_resids]},
                        style_lst=style_lst,
                        label_lst=[
                            [
                                f"Y{stick_resid}",
                                {
                                    "backgroundColor": "lightgray",
                                    "fontColor": "black",
                                    "backgroundOpacity": 0.5,
                                },
                                {"chain": chainid, "resi": stick_resid, "atom": "OH"},
                            ]
                        ],
                        cartoon_style=cartoon_style,
                        width=450,
                        height=300,
                        zoom=1.5,
                    )

                for col in table_col_lst:
                    table_df = (
                        pd.pivot_table(
                            data=rename_st_cols(
                                mask_equal(loop_df, loop_name, loop_conf)
                            ),
                            index=col,
                            columns=rename_col_dict[gene_class_col],
                            values=rename_col_dict[pdb_id_col],
                            aggfunc="nunique",
                            margins=True,
                        )
                        .reset_index()
                        .fillna("")
                    )

                    for col in list(table_df.columns):
                        table_df[col] = table_df[col].map(str)
                        table_df = fix_col(table_df, col)

                    show_st_table(table_df, st_col=loop_col)

        if nuc_class != gtp_name:
            st.markdown("---")

    write_st_end()