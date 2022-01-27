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

import numpy as np
from tqdm import tqdm
import itertools

from ..functions.table import order_rows
from ..functions.interf import calc_q_score
from ..functions.col import interf_cont_col, cb_dist_col
from ..functions.lst import str_to_lst
from ..functions.path import save_matrix


def calc_interf_dist(i, j, i_df, j_df):

    i_cont_lst = str_to_lst(i_df.at[i, interf_cont_col])
    j_cont_lst = str_to_lst(j_df.at[j, interf_cont_col])

    i_dist_lst = str_to_lst(i_df.at[i, cb_dist_col], return_float=True)
    j_dist_lst = str_to_lst(j_df.at[j, cb_dist_col], return_float=True)

    dist = calc_q_score(i_cont_lst, j_cont_lst, i_dist_lst, j_dist_lst)

    result = (i, j, dist)

    return result


def build_interf_matrix(fit_df, interf_matrix_path, pred_df=None):

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
        index_pairs, desc="Building interface matrix", position=0, leave=True
    ):

        result = calc_interf_dist(i_index, j_index, i_df, j_df)

        i_index = result[0]
        j_index = result[1]
        dist = result[2]

        matrix[i_index, j_index] = dist

        if pred_df is None:
            matrix[j_index, i_index] = dist

    save_matrix(interf_matrix_path, matrix)

    print("Built interface matrix!")