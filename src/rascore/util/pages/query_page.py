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

from ..functions import *
from ..scripts import *


def query_page():

    df = load_st_table(__file__)
    df = rename_st_cols(df)

    print(df.columns)

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

    for annot in annot_lst:
        mask_dict[rename_col_dict[annot]] = st.sidebar.multiselect(
            rename_col_dict[annot],
            ["All"] + lst_col(df, rename_col_dict[annot], unique=True),
        )

    st.markdown("#### Pivot Table")

    left_col, right_col = st.columns(2)

    mask_df = mask_st_table(df, mask_dict)

    row_lst = left_col.multiselect("Rows", [rename_col_dict[x] for x in annot_lst])

    col_lst = right_col.multiselect(
        "Columns",
        [rename_col_dict[x] for x in annot_lst if rename_col_dict[x] not in row_lst],
    )

    if len(row_lst) == 0 and len(col_lst) == 0:
        col_lst = rename_col_dict[sw1_col]
        row_lst = rename_col_dict[sw2_col]

    if len(row_lst) > 0 and len(col_lst) > 0:
        pivot_df = (
            pd.pivot_table(
                data=mask_df,
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

        download_st_df(pivot_df, "Download Pivot Table", "rascore_pivot.csv")

        st.markdown("---")

        st.markdown("#### Entries Table")

        show_st_dataframe(mask_df)

        download_st_df(mask_df, "Download Entries Table", "rascore_entries.csv")
