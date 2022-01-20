# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""
import numpy as np
from tqdm import tqdm
import scipy
import itertools

from ..functions import *


def calc_dist_dist(i_index, j_index, i_df, j_df):

    i_vect_1_col_lst = get_col_col_lst(i_df, vect_1_col)
    j_vect_1_col_lst = get_col_col_lst(j_df, vect_1_col)

    i_vect_2_col_lst = get_col_col_lst(i_df, vect_2_col)
    j_vect_2_col_lst = get_col_col_lst(j_df, vect_2_col)

    vect_1_col_lst = lst_inter(i_vect_1_col_lst, j_vect_1_col_lst)
    vect_2_col_lst = lst_inter(i_vect_2_col_lst, j_vect_2_col_lst)

    vect_col_lst = vect_1_col_lst + vect_2_col_lst

    vect_i_lst = tuple()
    vect_j_lst = tuple()

    for col in vect_col_lst:
        vect_i_lst += (str_to_lst(i_df.at[i_index, col]),)
        vect_j_lst += (str_to_lst(j_df.at[j_index, col]),)

    dist = scipy.spatial.distance.euclidean(vect_i_lst, vect_j_lst)

    return (i_index, j_index, dist)


def build_dist_matrix(
    fit_df,
    dist_matrix_path,
    pred_df=None,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

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
        index_pairs, desc="Building distance matrix", position=0, leave=True
    ):
        result = calc_dist_dist(
            i_index,
            j_index,
            i_df,
            j_df,
        )

        i_index = result[0]
        j_index = result[1]
        dist = result[2]

        matrix[i_index, j_index] = dist

        if pred_df is None:
            matrix[j_index, i_index] = dist

    save_matrix(dist_matrix_path, matrix)