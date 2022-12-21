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

from ..constants.conf import (
    loop_resid_dict,
    sw1_name,
    sw2_name,
    outlier_name,
    disorder_name,
    conf_color_dict,
)
from ..constants.nuc import nuc_color_dict, nuc_class_lst, nf_name, gdp_name, gtp_name

from ..functions.gui import (
    get_neighbor_path,
    load_st_table,
    get_html_text,
    write_st_end,
    show_st_structure,
    rename_st_cols,
    show_st_table,
    reorder_st_cols
)
from ..functions.path import (
    get_file_path,
    load_table,
    pages_str,
    data_str,
    rascore_str,
    cluster_str,
)
from ..functions.table import make_dict
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
    get_col_most_common,
    lst_col,
    fix_col,
)

reverse_col_dict = make_dict(list(rename_col_dict.values()),list(rename_col_dict.keys()))


def conformation_page():

    df = load_st_table(__file__)

    st.markdown("## Explore Conformations")
    st.markdown("---")

    st.sidebar.markdown(
        """
    **Note.** Only structures with completely modeled loops from our original clustering are displayed.
    Tables below structures include updated counts from *Rascore* database.
    """
    )

    data_path = get_neighbor_path(__file__, pages_str, data_str)

    sum_df = (
                    pd.pivot_table(
                        data=rename_st_cols(df),
                        index=rename_col_dict[sw2_name],
                        columns=rename_col_dict[sw1_name],
                        values=rename_col_dict[pdb_id_col],
                        aggfunc="nunique",
                        margins=True,
                    )
                    .fillna("")
                )

    for col in list(sum_df.columns):
        sum_df[col] = sum_df[col].map(str)
        sum_df = fix_col(sum_df, col)

    sum_df = reorder_st_cols(sum_df,sw2_name,sw1_name)
    sum_df = sum_df.reset_index()

    sum_df = sum_df.rename(columns={rename_col_dict[sw2_name]:'Conformation Label'})

    show_st_table(sum_df)

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

    nuc_name_dict = {nf_name: "Nucleotide-Free", gdp_name: "GDP-Bound", gtp_name: "GTP-Bound"}

    for nuc_class in nuc_class_lst:
        with st.expander(f"{nuc_class} Conformations ({nuc_name_dict[nuc_class]})", expanded=True):
            st.markdown(
                get_html_text(
                    {
                        nuc_class: nuc_color_dict[nuc_class],
                        f" Conformations ({nuc_name_dict[nuc_class]})": "#31333F",
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
                    loop_col = sw1_col

                elif loop_name == sw2_name:
                    cluster_df = sw2_df.copy(deep=True)
                    loop_col = sw2_col


                loop_df = mask_equal(nuc_df, loop_name, [x for x in lst_col(nuc_df, loop_name,unique=True) if outlier_name not in x and disorder_name not in x])

                conf_lst = get_col_most_common(loop_df, loop_name)

                if len(conf_lst) > 1:
                    conf_str = f"{loop_name} Conformations ({len(conf_lst)} in Total)"
                else:
                    conf_str = f"{loop_name} Conformation (Only 1)"

                loop_col.markdown(f"##### {conf_str}")

                loop_conf = loop_col.selectbox("Conformation Name", conf_lst)

                conf_df = mask_equal(loop_df,loop_name,loop_conf)

                pdb_id = loop_col.selectbox(
                    "PDB ID",
                    [
                        x.upper()
                        for x in lst_col(conf_df, pdb_id_col)
                        if x in lst_col(cluster_df, pdb_id_col)
                    ],
                )

                pdb_code = pdb_id[:4].lower()
                chainid = pdb_id[4:5]

                pdb_df = mask_equal(conf_df, pdb_id_col, f"{pdb_code}{chainid}")

                if loop_name == sw1_name:
                    sw1_color = conf_color_dict[loop_name][loop_conf]
                    sw2_color = "white"
                    stick_resids = 32
                elif loop_name == sw2_name:
                    sw1_color = "white"
                    sw2_color = conf_color_dict[loop_name][loop_conf]
                    stick_resids = 71

                show_st_structure(pdb_df,
                    stick_resids=stick_resids,
                    label_resids=True, 
                    zoom_resids=loop_resids,
                    sw1_color=sw1_color,
                    sw2_color=sw2_color,
                    zoom=1.5,
                    width=400,
                    height=300,
                    st_col=loop_col)

                for table_col in table_col_lst:
                    table_df = (
                        pd.pivot_table(
                            data=rename_st_cols(
                                mask_equal(loop_df, loop_name, loop_conf)
                            ),
                            index=table_col,
                            columns=rename_col_dict[gene_class_col],
                            values=rename_col_dict[pdb_id_col],
                            aggfunc="nunique",
                            margins=True,
                        )
                        .fillna("")
                    )

                    for col in list(table_df.columns):
                        table_df[col] = table_df[col].map(str)
                        table_df = fix_col(table_df, col)

                    table_df = reorder_st_cols(table_df,reverse_col_dict[table_col],gene_class_col)
                    table_df = table_df.reset_index()

                    show_st_table(table_df, st_col=loop_col)

    write_st_end()