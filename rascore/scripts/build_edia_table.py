# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import pandas as pd
from tqdm import tqdm

from functions import *


def get_index_edia(
    df,
    index,
    edia_dict,
    resid_lst,
    max_ca_dist=2,
    ext_mult=1,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    index_df = pd.DataFrame()

    pdb_code = df.at[index, pdb_code_col]

    if pdb_code in list(edia_dict.keys()):
        chainid = df.at[index, chainid_col]
        if chainid in list(edia_dict[pdb_code].keys()):
            coord_path = df.at[index, coord_path_col]
            modelid = int(df.at[index, modelid_col])

            edia_resid_lst = list(edia_dict[pdb_code][chainid].keys())

            add_resid_lst = build_add_resid_lst(
                coord_path,
                modelid,
                chainid,
                resid_lst,
                edia_resid_lst,
                max_ca_dist=max_ca_dist,
                ext_mult=ext_mult,
            )

            i = 0

            for n, add_resid in enumerate(add_resid_lst):

                if str(add_resid) in edia_resid_lst:

                    for atomid in list(
                        edia_dict[pdb_code][chainid][str(add_resid)].keys()
                    ):

                        index_df.at[i, pdb_code_col] = pdb_code
                        index_df.at[i, chainid_col] = chainid
                        index_df.at[i, resid_col] = n + 1
                        index_df.at[i, atomid_col] = atomid
                        index_df.at[i, edia_col] = edia_dict[pdb_code][chainid][
                            str(add_resid)
                        ][atomid][edia_col]
                        index_df.at[i, b_factor_col] = edia_dict[pdb_code][chainid][
                            str(add_resid)
                        ][atomid][b_factor_col]

                        i += 1

    return index_df


def build_edia_table(
    df, edia_dict, edia_resids, edia_table_path=None, coord_path_col=None
):

    resid_lst = res_to_lst(edia_resids)

    edia_df = pd.DataFrame()

    for index in tqdm(
        list(df.index.values), desc="Building EDIA table", position=0, leave=True
    ):

        edia_df = pd.concat(
            [
                edia_df,
                get_index_edia(
                    df, index, edia_dict, resid_lst, coord_path_col=coord_path_col
                ),
            ],
            sort=False,
        )

    edia_df = merge_tables(left_df=edia_df, right_df=df)

    edia_df = edia_df.reset_index(drop=True)

    if edia_table_path is not None:
        save_table(edia_table_path, edia_df)
    else:
        return edia_df

    print("Built EDIA table!")
