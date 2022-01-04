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


def get_index_pocket(
    df,
    index,
    pocket_dict,
    min_volume=None,
    min_score=None,
    search_cont_lst=None,
    search_max_dist=0.7,
    coord_path_col=None,
    use_simpson=False,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    index_df = get_df_at_index(df, index)

    coord_path = index_df.at[index, coord_path_col]
    coord_path = coord_path.replace(".cif", ".pdb")

    pocket_df = pd.DataFrame()

    for pocket in list(pocket_dict[coord_path].keys()):

        add_pocket = True

        if min_volume is not None:
            if float(pocket_dict[coord_path][pocket][pocket_volume_col]) < min_volume:
                add_pocket = False

        if min_score is not None:
            if float(pocket_dict[coord_path][pocket][pocket_score_col]) < min_score:
                add_pocket = False

        if search_cont_lst is not None:
            pocket_cont_lst = str_to_lst(
                pocket_dict[coord_path][pocket][pocket_cont_col]
            )
            if use_simpson:
                pocket_dist = calc_simpson(
                    pocket_cont_lst,
                    search_cont_lst,
                    return_dist=True,
                )

            else:
                pocket_dist = calc_jaccard(
                    pocket_cont_lst,
                    search_cont_lst,
                    return_dist=True,
                )

            if pocket_dist > search_max_dist:
                add_pocket = False

        if add_pocket:
            temp_df = index_df.copy(deep=True)

            for col in list(pocket_dict[coord_path][pocket].keys()):
                temp_df.at[index, col] = pocket_dict[coord_path][pocket][col]

            temp_df.at[index, pocket_col] = pocket
            temp_df.at[index, pocket_lig_col] = "STP"

            pocket_df = pd.concat([pocket_df, temp_df], sort=False)

    return pocket_df


def build_pocket_table(
    df,
    pocket_dict,
    pocket_table_path=None,
    min_volume=None,
    min_score=None,
    search_cont_lst=None,
    search_max_dist=0.7,
    use_simpson=False,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    pocket_df = pd.DataFrame()

    for index in tqdm(
        list(df.index.values), desc="Building pocket table", position=0, leave=True
    ):

        pocket_df = pd.concat(
            [
                pocket_df,
                get_index_pocket(
                    df,
                    index,
                    pocket_dict,
                    min_volume=min_volume,
                    min_score=min_score,
                    search_cont_lst=search_cont_lst,
                    search_max_dist=search_max_dist,
                    use_simpson=use_simpson,
                    coord_path_col=coord_path_col,
                ),
            ],
            sort=False,
        )

    pocket_df = pocket_df.reset_index(drop=True)

    if pocket_table_path is not None:
        save_table(pocket_table_path, pocket_df)
    else:
        return pocket_df

    print("Built pocket table!")
