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
import numpy as np
import math
import itertools

from .lst import move_end_lst
from .col import (
    cluster_col,
    pdb_code_col,
    cf_col,
    nn_dist_col,
    constr_dist_col,
    entropy_col,
    common_col,
    occup_col,
    rama_col,
    rotamer_col,
    bb_seq_col,
    sc_seq_col,
    iso_col,
    dist_col_lst,
    interf_area_col,
    pocket_volume_col,
    pocket_score_col,
    mean_col,
    max_col,
)
from .table import (
    build_count_table,
    merge_tables,
    mask_equal,
    build_count_dict,
    build_col_count_dict,
    lst_col,
    lst_by_freq,
    get_val_index_lst,
    mask_matrix,
    get_col_most_common,
)


def is_noise(val):

    noise_status = False
    if type(val) == str:
        if 'Noise' in val:
            noise_status = True
    elif val == -1:
        noise_status = True

    return noise_status


def label_cluster(val, cluster_order):

    if is_noise(val):
        return "Noise"
    else:
        return str(cluster_order.index(val) + 1)


def order_clusters(df, order_col=None):

    temp_df = df.copy(deep=True)
    if order_col is not None:
        temp_df = temp_df.loc[:, [cluster_col, order_col]]
        temp_df = temp_df.drop_duplicates()

    cluster_lst = lst_col(temp_df, cluster_col)

    cluster_order = lst_by_freq(cluster_lst)

    for cluster in cluster_order:
        if is_noise(cluster):
            cluster_order = move_end_lst(cluster_order, "Noise")

    df[cluster_col] = df.apply(
        lambda x: label_cluster(x[cluster_col], cluster_order),
        axis=1,
    )

    return df


def renumber_clusters(df, min_samples, min_pdb=2, min_cf=1, order_col=None):

    count_dict = build_count_dict(lst_col(df, cluster_col, return_str=True))

    check_pdb = False
    if pdb_code_col in list(df.columns):
        pdb_dict = build_col_count_dict(
            df, cluster_col, col_lst=[pdb_code_col], return_str=True
        )
        check_pdb = True

    check_cf = False
    if cf_col in list(df.columns):
        cf_dict = build_col_count_dict(
            df, cluster_col, col_lst=[cf_col], return_str=True
        )
        check_cf = True

    df[cluster_col] = df[cluster_col].map(str)

    label_dict = dict()
    for cluster in lst_col(df, cluster_col, unique=True, return_str=True):

        label = cluster

        if count_dict[cluster] < min_samples:
            label = "Noise"

        if check_pdb:
            if pdb_dict[cluster] < min_pdb:
                label = "Noise"

        if check_cf:
            if cf_dict[cluster] < min_cf:
                label = "Noise"

        label_dict[cluster] = label

    df[cluster_col] = df[cluster_col].map(label_dict)

    df = order_clusters(df, order_col=order_col)

    return df


def get_index_dist_lst(matrix, index_lst, index=None):

    if index is not None:
        index_lst.remove(index)

    return matrix[index_lst, index]


def calc_dist_stat(dist_lst, method="mean", dec=2):

    if len(dist_lst) < 2:
        dist = dist_lst[0]
    else:
        if method == "mean":
            dist = np.mean(dist_lst)
        elif method == "max":
            dist = np.max(dist_lst)
        elif method == "min":
            dist = np.min(dist_lst)
        elif method == "median":
            dist = np.median(dist_lst)
        elif method == "std":
            dist = np.std(dist_lst)

    return round(dist, dec)


def calc_cutoff(dist_lst, dec=2):

    return round(
        calc_dist_stat(dist_lst, method="mean")
        + (2 * calc_dist_stat(dist_lst, method="std")),
        dec,
    )


def build_cutoff_dict(df, matrix, method="mean"):

    cluster_lst = lst_col(df, cluster_col, unique=True, return_str=True)

    cutoff_dict = dict()

    for cluster in cluster_lst:

        index_lst = get_val_index_lst(df, cluster_col, cluster)

        dist_lst = list()

        if len(index_lst) == 1:
            dist_lst.append(0.0)
        else:
            for index in index_lst:
                dist_lst.append(
                    calc_dist_stat(
                        get_index_dist_lst(matrix, index_lst, index=index),
                        method=method,
                    )
                )

        cutoff_dict[cluster] = calc_cutoff(dist_lst)

    return cutoff_dict


def apply_dist_constr(
    df,
    dist_col,
    matrix,
    min_samples=None,
    max_dist=None,
    constr_method=None,
    order_col=None,
):
    df[cluster_col] = df[cluster_col].map(str)

    df_index_lst = list(df.index.values)

    if max_dist is None:
        cutoff_dict = build_cutoff_dict(df, matrix, method=constr_method)

    pre_count_dict = build_col_count_dict(df, cluster_col)

    calc_dist = True
    while calc_dist:

        calc_dist = False

        for index in df_index_lst:

            cluster = df.at[index, cluster_col]

            index_lst = get_val_index_lst(df, cluster_col, cluster)

            if min_samples is not None:
                if not is_noise(cluster):
                    if len(index_lst) < min_samples:
                        cluster = "Noise"
                        df.at[index, cluster_col] = cluster
                        calc_dist = True

            if len(index_lst) == 1:
                dist = 0.0
            else:
                dist = calc_dist_stat(
                    get_index_dist_lst(matrix, index_lst, index=index),
                    method=constr_method,
                )

            df.at[index, dist_col] = dist

            if not is_noise(cluster):

                prune_cluster = False

                dist_cutoff = max_dist
                if max_dist is None:
                    dist_cutoff = cutoff_dict[cluster]

                if dist > dist_cutoff:
                    prune_cluster = True

                if prune_cluster:
                    df.at[index, cluster_col] = "Noise"
                    calc_dist = True

    post_count_dict = build_col_count_dict(df, cluster_col)

    pruned = 0
    for cluster in list(post_count_dict.keys()):
        if cluster != "Noise":
            pruned += pre_count_dict[cluster] - post_count_dict[cluster]

    df = order_clusters(df, order_col=order_col)

    return df, pruned


def prune_cluster_members(
    df,
    min_samples,
    matrix,
    max_nn_dist=None,
    constr_matrix=None,
    max_constr_dist=None,
    constr_method="mean",
    min_pdb=3,
    min_cf=1,
    order_col=None,
):

    df = renumber_clusters(
        df, min_samples, min_pdb=min_pdb, min_cf=min_cf, order_col=order_col
    )

    df, nn_pruned = apply_dist_constr(
        df,
        nn_dist_col,
        matrix,
        min_samples=min_samples,
        max_dist=max_nn_dist,
        constr_method="min",
        order_col=order_col,
    )

    if constr_matrix is not None:
        if max_constr_dist is None:
            max_constr_dist = build_cutoff_dict(df, constr_matrix, method=constr_method)

        df, constr_pruned = apply_dist_constr(
            df,
            constr_dist_col,
            constr_matrix,
            min_samples=min_samples,
            max_dist=max_constr_dist,
            constr_method=constr_method,
            order_col=order_col,
        )

    df = renumber_clusters(
        df,
        min_samples,
        min_pdb=min_pdb,
        min_cf=min_cf,
        order_col=order_col,
    )

    if constr_matrix is not None:
        return (df, nn_pruned, constr_pruned)
    else:
        return (df, nn_pruned)


def merge_clusters(df, matrix, max_merge_dist, merge_method="mean"):

    cluster_lst = [x for x in lst_col(df, cluster_col, unique=True) if not is_noise(x)]

    cluster_pairs = itertools.combinations(cluster_lst, 2)

    merge_dict = dict()

    for cluster_1, cluster_2 in cluster_pairs:

        index_1_lst = get_val_index_lst(df, cluster_col, cluster_1)
        index_2_lst = get_val_index_lst(df, cluster_col, cluster_2)

        if (
            calc_dist_stat(
                mask_matrix(matrix, index_1_lst, index_2_lst), method=merge_method
            )
            <= max_merge_dist
        ):
            merge_cluster = cluster_lst[
                min(cluster_lst.index(cluster_1), cluster_lst.index(cluster_2))
            ]
            merge_dict[cluster_1] = merge_cluster
            merge_dict[cluster_2] = merge_cluster

    df[cluster_col] = df[cluster_col].replace(merge_dict)

    return df


def calc_entropy(vals, dec=2):

    return round(
        (
            float(-1)
            * float(np.max(vals))
            / float(np.sum(vals))
            * (math.log10(float(np.max(vals)) / float(np.sum(vals))))
        ),
        dec,
    )


def calc_occupancy(vals, dec=2):

    return round(float(np.max(vals)) / float(np.sum(vals)), dec)


def add_dih_stats(temp_df, sum_df, col, index):

    count_dict = build_col_count_dict(temp_df, col)

    common = [x for x in get_col_most_common(temp_df, col) if "-" not in x]
    if len(common) == 0:
        common = "-"
    else:
        common = common[0]
    entropy = calc_entropy(list(count_dict.values()))
    occupancy = calc_occupancy(list(count_dict.values()))

    sum_df.at[index, f"{common_col}_{col}"] = common
    sum_df.at[index, f"{entropy_col}_{col}"] = entropy
    sum_df.at[index, f"{occup_col}_{col}"] = occupancy

    return sum_df


def dist_to_dih(dist, dec=2):

    return round(math.degrees(math.acos(1 - (float(dist) / 2))), dec)


def build_sum_table(df):

    df_col_lst = list(df.columns)

    dih_dist = False
    if rama_col in df_col_lst or rotamer_col in df_col_lst:
        dih_dist = True

    add_col_lst = [rama_col, rotamer_col, bb_seq_col, sc_seq_col, iso_col]

    df[cluster_col] = df[cluster_col].map(str)

    cluster_lst = lst_col(df, cluster_col, unique=True, return_str=True)

    sum_df = pd.DataFrame()

    stat_col_lst = (
        dist_col_lst + [interf_area_col] + [pocket_volume_col] + [pocket_score_col]
    )

    for cluster in cluster_lst:

        temp_df = mask_equal(df, cluster_col, cluster)

        index = cluster_lst.index(cluster)

        sum_df.at[index, cluster_col] = cluster

        if "Noise" not in cluster and "Outlier" not in cluster and "Disordered" not in cluster:
            for col in stat_col_lst:
                if col in df_col_lst:
                    dist_lst = lst_col(temp_df, col, return_float=True)

                    mean_dist = calc_dist_stat(dist_lst, method="mean")
                    max_dist = calc_dist_stat(dist_lst, method="max")

                    val_lst = [mean_dist, max_dist]

                    if dih_dist:
                        if col == nn_dist_col:
                            val_lst = [dist_to_dih(x) for x in val_lst]

                    sum_df.at[index, f"{col}_{mean_col}"] = val_lst[0]
                    sum_df.at[index, f"{col}_{max_col}"] = val_lst[1]

            for add_col in add_col_lst:
                if add_col in df_col_lst:
                    sum_df = add_dih_stats(temp_df, sum_df, add_col, index)

    sum_df = merge_tables(sum_df, build_count_table(df, cluster_col))

    sum_df = sum_df.fillna("-")

    return sum_df