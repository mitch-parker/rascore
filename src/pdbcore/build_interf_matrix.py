# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import numpy as np
from tqdm import tqdm
import itertools
import math

from util.data import order_rows
from util.lst import str_to_lst, lst_unique
from util.path import save_matrix
from util.col import interf_cont_col, cb_dist_col


def calc_interf_dist(i, j, i_df, j_df):

    i_cont_lst = str_to_lst(i_df.at[i, interf_cont_col])
    j_cont_lst = str_to_lst(j_df.at[j, interf_cont_col])

    i_dist_lst = str_to_lst(i_df.at[i, cb_dist_col], return_float=True)
    j_dist_lst = str_to_lst(j_df.at[j, cb_dist_col], return_float=True)

    cont_lst = lst_unique(i_cont_lst, j_cont_lst)

    q_sum = 0
    w_sum = 0

    for cont in cont_lst:

        i_cb_dist = 12
        if cont in i_cont_lst:
            i_cb_dist = i_dist_lst[i_cont_lst.index(cont)]

        j_cb_dist = 12
        if cont in j_cont_lst:
            j_cb_dist = j_dist_lst[j_cont_lst.index(cont)]

        min_cb_diff = min(i_cb_dist, j_cb_dist)
        abs_cb_diff = abs(i_cb_dist - j_cb_dist)

        if min_cb_diff > 5:
            weight = math.exp(-((min_cb_diff - 5) ** 2) / 9.158)
        else:
            weight = 1

        q_sum += weight * math.exp(-0.5 * abs_cb_diff)
        w_sum += weight

    if w_sum == 0:
        dist = 1
    else:
        dist = 1 - (q_sum / w_sum)

    result = (i, j, dist)

    return result


def build_interf_matrix(included_df, interf_matrix_path, removed_df=None):

    included_df = order_rows(included_df)

    included_df = order_rows(included_df)
    j_df = included_df.copy(deep=True)

    if removed_df is None:
        i_df = included_df.copy(deep=True)
    else:
        removed_df = order_rows(removed_df)
        i_df = removed_df.copy(deep=True)

    i_index_lst = list(i_df.index.values)
    j_index_lst = list(j_df.index.values)

    if removed_df is not None:
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

        if removed_df is None:
            matrix[j_index, i_index] = dist

    save_matrix(interf_matrix_path, matrix)

    print("Built interface matrix!")