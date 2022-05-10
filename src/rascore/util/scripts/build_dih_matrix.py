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
import math

from ..functions.lst import type_lst, lst_inter
from ..functions.table import extract_int, get_col_val_lst, get_val_col, order_rows
from ..functions.col import phi_col, psi_col, dih_col_lst
from ..functions.path import save_matrix


def calc_flip_dist(norm_dist, resid_dict, index_dict, norm_dict, flip_dict, mean=False):

    if norm_dist <= 2:
        flip_dist = norm_dist
    else:
        flip_dist_lst = list()
        bb_resid_lst = resid_dict[phi_col]
        min_resid = min(bb_resid_lst)
        max_resid = max(bb_resid_lst)
        for flip_diff in list(flip_dict.keys()):
            for index, resid in enumerate(bb_resid_lst):
                if resid != min_resid and resid != max_resid:
                    curr_resid = resid
                    next_resid = bb_resid_lst[index + 1]
                    if abs(extract_int(next_resid) - extract_int(curr_resid)) < 2:
                        if curr_resid in list(
                            index_dict[psi_col].keys()
                        ) and next_resid in list(index_dict[phi_col].keys()):
                            curr_index = index_dict[psi_col][curr_resid]
                            next_index = index_dict[phi_col][next_resid]

                            temp_psi_vals = np.copy(norm_dict[psi_col])
                            temp_phi_vals = np.copy(norm_dict[phi_col])

                            temp_psi_vals[curr_index] = flip_dict[flip_diff][psi_col][
                                curr_index
                            ]
                            temp_phi_vals[next_index] = flip_dict[flip_diff][phi_col][
                                next_index
                            ]

                            temp_flip_vals = np.array([])

                            temp_flip_vals = np.concatenate(
                                (temp_flip_vals, temp_psi_vals), axis=0
                            )
                            temp_flip_vals = np.concatenate(
                                (temp_flip_vals, temp_phi_vals), axis=0
                            )

                            for dih_col in dih_col_lst:
                                if dih_col != psi_col and dih_col != phi_col:
                                    norm_angle_vals = norm_dict[dih_col]
                                    if len(norm_angle_vals) > 0:
                                        temp_flip_vals = np.concatenate(
                                            (temp_flip_vals, norm_angle_vals), axis=0
                                        )
                            if mean:
                                temp_flip_dist = np.mean(temp_flip_vals)
                            else:
                                temp_flip_dist = np.max(temp_flip_vals)
                            flip_dist_lst.append(temp_flip_dist)

        if len(flip_dist_lst) == 0:
            flip_dist = norm_dist
        else:
            flip_dist = min(min(flip_dist_lst), norm_dist)

    return flip_dist


def calc_dih_dist(i, j, i_df, j_df, flip_diff=180):

    resid_dict = dict()

    for dih_col in dih_col_lst:

        resid_dict[dih_col] = lst_inter(get_col_val_lst(i_df, dih_col), get_col_val_lst(j_df, dih_col))

    norm_dict = dict()
    flip_dict = dict()

    for dih_col in dih_col_lst:

        norm_dict[dih_col] = np.array([])

    flip_diff_lst = type_lst(flip_diff)

    for flip_diff in flip_diff_lst:
        flip_dict[flip_diff] = {phi_col: np.array([]), psi_col: np.array([])}

    index_dict = dict()

    for dih_col in dih_col_lst:

        resid_lst = resid_dict[dih_col]

        index_dict[dih_col] = dict()
        index = 0

        for resid in resid_lst:

            val_col = get_val_col(dih_col, resid)

            norm_val_i = float(i_df.at[i, val_col])
            norm_val_j = float(j_df.at[j, val_col])

            if norm_val_i != 999.00 and norm_val_j != 999.00:

                norm_dist = 2 * (
                    1 - math.cos((math.radians(norm_val_i) - math.radians(norm_val_j)))
                )

                norm_dict[dih_col] = np.append(norm_dict[dih_col], norm_dist)

                for flip_diff in flip_diff_lst:

                    if dih_col == psi_col or dih_col == phi_col:
                        flip_dist = 2 * (
                            1
                            - math.cos(
                                (
                                    math.radians(norm_val_i + flip_diff)
                                    - math.radians(norm_val_j)
                                )
                            )
                        )

                        flip_dict[flip_diff][dih_col] = np.append(
                            flip_dict[flip_diff][dih_col], flip_dist
                        )

                index_dict[dih_col][resid] = index
                index += 1

    norm_vals = np.array([])

    for dih_col in dih_col_lst:
        norm_vals = np.concatenate((norm_vals, norm_dict[dih_col]), axis=0)

    if len(norm_vals) == 0:
        max_norm_dist = 4.0
        mean_norm_dist = 4.0
        max_flip_dist = 4.0
        mean_flip_dist = 4.0
    else:
        max_norm_dist = np.max(norm_vals)
        mean_norm_dist = np.mean(norm_vals)
        max_flip_dist = calc_flip_dist(
            max_norm_dist, resid_dict, index_dict, norm_dict, flip_dict, mean=False
        )
        mean_flip_dist = calc_flip_dist(
            mean_norm_dist, resid_dict, index_dict, norm_dict, flip_dict, mean=True
        )

    max_norm_result = (i, j, max_norm_dist)
    mean_norm_result = (i, j, mean_norm_dist)
    max_flip_result = (i, j, max_flip_dist)
    mean_flip_result = (i, j, mean_flip_dist)

    return (max_norm_result, mean_norm_result, max_flip_result, mean_flip_result)


def build_dih_matrix(
    fit_df,
    max_norm_path,
    mean_norm_path=None,
    max_flip_path=None,
    mean_flip_path=None,
    pred_df=None,
    flip_diff=180,
):

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

    max_norm_matrix = np.zeros((len(i_index_lst), len(j_index_lst)))
    mean_norm_matrix = np.zeros((len(i_index_lst), len(j_index_lst)))
    max_flip_matrix = np.zeros((len(i_index_lst), len(j_index_lst)))
    mean_flip_matrix = np.zeros((len(i_index_lst), len(j_index_lst)))

    matrix_lst = [max_norm_matrix, mean_norm_matrix, max_flip_matrix, mean_flip_matrix]

    for i_index, j_index in tqdm(
        index_pairs, desc="Building dihedral matrix", position=0, leave=True
    ):
        dih_result = calc_dih_dist(
            i_index,
            j_index,
            i_df,
            j_df,
            flip_diff=flip_diff,
        )

        result_lst = [
            dih_result[0],
            dih_result[1],
            dih_result[2],
            dih_result[3],
        ]

        for result, matrix in zip(result_lst, matrix_lst):

            i_index = result[0]
            j_index = result[1]
            dist = result[2]

            matrix[i_index, j_index] = dist

            if pred_df is None:
                matrix[j_index, i_index] = dist

    save_matrix(max_norm_path, max_norm_matrix)

    if mean_norm_path is not None:
        save_matrix(mean_norm_path, mean_norm_matrix)

    if max_flip_path is not None:
        save_matrix(max_flip_path, max_flip_matrix)

    if mean_flip_path is not None:
        save_matrix(mean_flip_path, mean_flip_matrix)

    print("Built dihedral matrix!")
