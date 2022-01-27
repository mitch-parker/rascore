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
import numpy as np
from tqdm import tqdm

from ..functions.col import (
    core_path_col,
    pocket_volume_col,
    pocket_score_col,
    pocket_cont_col,
    pocket_col,
    pocket_id_col,
    pdb_id_col,
)
from ..functions.table import get_df_at_index, fix_col
from ..functions.lst import str_to_lst, calc_jaccard, calc_simpson
from ..functions.path import save_table


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

    if search_cont_lst is not None:
        if not any(isinstance(x, list) for x in search_cont_lst):
            search_cont_lst = [search_cont_lst]
            for i, search_cont in enumerate(search_cont_lst):
                search_cont_lst[i] = [str(x) for x in search_cont]

    index_df = get_df_at_index(df, index)

    coord_path = index_df.at[index, coord_path_col]
    coord_path = coord_path.replace(".cif", ".pdb")

    pocket_df = pd.DataFrame()

    if coord_path in list(pocket_dict.keys()):
        for pocket in list(pocket_dict[coord_path].keys()):

            add_pocket = True

            if min_volume is not None:
                if (
                    float(pocket_dict[coord_path][pocket][pocket_volume_col])
                    < min_volume
                ):
                    add_pocket = False

            if min_score is not None:
                if float(pocket_dict[coord_path][pocket][pocket_score_col]) < min_score:
                    add_pocket = False

            if search_cont_lst is not None:

                pocket_cont = str_to_lst(
                    pocket_dict[coord_path][pocket][pocket_cont_col]
                )

                pocket_dist_lst = list()
                for search_cont in search_cont_lst:
                    if use_simpson:
                        pocket_dist = calc_simpson(
                            pocket_cont,
                            search_cont,
                            return_dist=True,
                        )

                    else:
                        pocket_dist = calc_jaccard(
                            pocket_cont,
                            search_cont,
                            return_dist=True,
                        )
                    pocket_dist_lst.append(pocket_dist)

                if np.mean(pocket_dist) > search_max_dist:
                    add_pocket = False

            if add_pocket:
                temp_df = index_df.copy(deep=True)

                for col in list(pocket_dict[coord_path][pocket].keys()):
                    temp_df.at[index, col] = pocket_dict[coord_path][pocket][col]

                temp_df.at[index, pocket_col] = pocket

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

    pocket_df[pocket_col] = pocket_df[pocket_col].map(str)
    pocket_df = fix_col(pocket_df, pocket_col)

    pocket_df[pocket_id_col] = pocket_df[pdb_id_col] + pocket_df[pocket_col]

    if pocket_table_path is not None:
        save_table(pocket_table_path, pocket_df)
    else:
        return pocket_df

    print("Built pocket table!")
