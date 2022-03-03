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

import numpy as np
from tqdm import tqdm
import itertools

from ..functions.table import order_rows
from ..functions.lst import str_to_lst, calc_simpson, calc_jaccard
from ..functions.col import pocket_cont_col
from ..functions.path import save_matrix


def build_pocket_matrix(fit_df, pocket_matrix_path, pred_df=None, use_simpson=False):

    fit_df = order_rows(fit_df)

    fit_df = order_rows(fit_df)
    j_df = fit_df.copy(deep=True)

    if pred_df is None:
        i_df = fit_df.copy(deep=True)
    else:
        pred_df = order_rows(pred_df)
        i_df = pred_df.copy(deep=True)

    i_index_lst = list(i_df.index.values)
    j_index_lst = list(j_df.index.values)

    if pred_df is not None:
        index_pairs = itertools.product(i_index_lst, j_index_lst)
    else:
        index_pairs = itertools.combinations(i_index_lst, 2)

    matrix = np.zeros((len(i_index_lst), len(j_index_lst)))

    for i_index, j_index in tqdm(
        index_pairs, desc="Building pocket matrix", position=0, leave=True
    ):

        i_cont_lst = str_to_lst(i_df.at[i_index, pocket_cont_col])
        j_cont_lst = str_to_lst(j_df.at[j_index, pocket_cont_col])

        if use_simpson:
            dist = calc_simpson(
                i_cont_lst,
                j_cont_lst,
                return_dist=True,
            )

        else:
            dist = calc_jaccard(
                i_cont_lst,
                j_cont_lst,
                return_dist=True,
            )

        matrix[i_index, j_index] = dist

        if pred_df is None:
            matrix[j_index, i_index] = dist

    save_matrix(pocket_matrix_path, matrix)

    print("Built pocket matrix!")