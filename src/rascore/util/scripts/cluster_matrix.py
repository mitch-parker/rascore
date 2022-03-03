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
from tqdm import tqdm
import itertools
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples

from ..functions.table import mask_matrix, mask_unequal, build_col_count_dict, lst_col
from ..functions.lst import (
    calc_simpson,
    lst_nums,
    build_range_lst,
    sort_lst,
    lst_unique,
    get_lst_val_indices,
)
from ..functions.col import (
    cluster_col,
    silh_col,
    simi_col,
    total_complete_col,
    cluster_count_col,
    total_cluster_pre_col,
    total_cluster_post_col,
    total_noise_pre_col,
    total_noise_post_col,
    total_pruned_nn_col,
    total_pruned_constr_col,
    select_col,
)
from ..functions.cluster import (
    merge_clusters,
    build_cutoff_dict,
    prune_cluster_members,
    build_sum_table,
)
from ..functions.path import save_table


def run_dbscan(matrix, eps, min_samples):

    return (
        DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
        .fit(matrix)
        .labels_
    )


def build_pass_lst(label_lst, matrix, max_dist, min_silh):

    pass_lst = list()

    cluster_lst = list(set(label_lst))

    if len(cluster_lst) > 1:

        silhouette_lst = silhouette_samples(matrix, label_lst, metric="precomputed")

        for cluster in cluster_lst:
            if cluster != -1:

                cluster_index_lst = [
                    index for index, label in enumerate(label_lst) if label == cluster
                ]

                max_cluster_dist = np.max(
                    mask_matrix(matrix, cluster_index_lst, cluster_index_lst)
                )
                mean_cluster_silh = np.mean(silhouette_lst[cluster_index_lst])

                if max_cluster_dist <= max_dist and mean_cluster_silh > min_silh:
                    pass_lst.append(tuple([label_lst, cluster_index_lst]))

    return pass_lst


def check_connect(compare_lst, i_index, j_index, min_simi):

    overlap = calc_simpson(compare_lst[i_index][1], compare_lst[j_index][1])

    if overlap >= min_simi:
        status = 1

    if overlap < min_simi:
        status = 0

    return tuple([i_index, j_index, status])


def connect_graph(graph_dict, node_lst, result_dict):

    for node in node_lst:

        subgraph = set(
            [k for k, v in list(graph_dict.items()) if graph_dict[node].intersection(v)]
        )

        if subgraph:

            node_lst = [x for x in node_lst if x not in list(subgraph)]
            node_lst = node_lst
            result_dict[node] = subgraph

            result_dict = connect_graph(graph_dict, node_lst, result_dict)

        else:
            continue

        if len(node_lst) == 1 or len(node_lst) == 0:
            return result_dict
        break

    return result_dict


def run_grid_graph(df, result_lst, matrix, min_samples, min_silh, min_simi, max_dist):

    compare_lst = list()

    for label_lst in result_lst:

        pass_lst = build_pass_lst(label_lst, matrix, max_dist, min_silh)

        compare_lst += pass_lst

    total_comparisons = len(compare_lst)

    index_lst = lst_nums(0, total_comparisons - 1)

    connect_matrix = np.zeros((total_comparisons, total_comparisons))

    index_pairs = itertools.combinations(index_lst, 2)

    for i_index, j_index in index_pairs:

        result = check_connect(compare_lst, i_index, j_index, min_simi)

        i_index = result[0]
        j_index = result[1]
        status = result[2]

        connect_matrix[i_index, j_index] = status
        connect_matrix[j_index, i_index] = status

    node_lst = list(range(len(connect_matrix)))

    graph_dict = dict()
    result_dict = dict()

    for node in node_lst:

        neighbors = list(np.where(connect_matrix[:][node] == 1)[0])
        edges = neighbors + [node]
        graph_dict[node] = set(edges)

    result_dict = connect_graph(graph_dict, node_lst, result_dict)

    cluster_dict = dict()
    cluster_label = 0

    for node_lst in list(result_dict.values()):

        final_index_lst = set()

        for node in node_lst:

            index_lst = set(compare_lst[node][1])
            final_index_lst = final_index_lst.union(index_lst)

        if len(final_index_lst) >= min_samples:
            cluster_dict[cluster_label] = final_index_lst
            cluster_label = cluster_label + 1

    for cluster_label in list(cluster_dict.keys()):
        for index in cluster_dict[cluster_label]:
            df.at[index, cluster_col] = cluster_label

    if cluster_col in list(df.columns):
        df[cluster_col] = df[cluster_col].fillna(-1)
    else:
        df[cluster_col] = -1

    return df


def cluster_matrix(
    df,
    matrix,
    cluster_table_path,
    sum_table_path=None,
    report_table_path=None,
    eps_range="0.1-1.6",
    eps_step=0.1,
    min_samples_range="5-15",
    min_samples_step=1,
    silh_range=0.6,
    silh_step=0.1,
    simi_range=0.9,
    simi_step=0.1,
    max_dist=3.75,
    max_nn_dist=None,
    min_pdb=3,
    min_cf=1,
    min_min_samples=None,
    constr_matrix=None,
    max_constr_dist=None,
    merge_constr_dist=None,
    constr_method="mean",
    cluster_lim=None,
):

    eps_lst = build_range_lst(eps_range, eps_step, type=float)
    min_samples_lst = build_range_lst(min_samples_range, min_samples_step, type=int)
    silh_lst = build_range_lst(silh_range, silh_step, type=float)
    simi_lst = build_range_lst(simi_range, simi_step, type=float)

    if min_min_samples is None:
        min_min_samples = min(min_samples_lst)

    result_lst = list()
    for eps in tqdm(eps_lst, desc="Running Clustering", position=0, leave=True):
        for min_samples in min_samples_lst:
            result_lst.append(run_dbscan(matrix, eps, min_samples))

    run_lst = list()
    for min_silh in silh_lst:
        for min_simi in simi_lst:
            run_lst.append((min_silh, min_simi))

    if report_table_path is not None:
        report_df = pd.DataFrame()
        i = 0

    grid_graph_lst = list()
    cluster_count_lst = list()
    total_noise_lst = list()

    for run in tqdm(run_lst, desc="Running GridGraph", position=0, leave=True):

        min_silh = run[0]
        min_simi = run[1]

        df = run_grid_graph(
            df,
            result_lst,
            matrix,
            min_min_samples,
            min_silh,
            min_simi,
            max_dist,
        )

        if constr_matrix is not None:
            if merge_constr_dist is not None:
                if max_constr_dist is None:
                    max_constr_dist = build_cutoff_dict(
                        df, constr_matrix, method=constr_method
                    )
                df = merge_clusters(
                    df, constr_matrix, merge_constr_dist, merge_method=constr_method
                )

        result = prune_cluster_members(
            df,
            min_min_samples,
            matrix,
            max_nn_dist=max_nn_dist,
            constr_matrix=constr_matrix,
            max_constr_dist=max_constr_dist,
            constr_method=constr_method,
            min_pdb=min_pdb,
            min_cf=min_cf,
        )

        df = result[0]
        nn_pruned = result[1]
        total_pruned = nn_pruned
        if constr_matrix is not None:
            constr_pruned = result[2]
            total_pruned += constr_pruned

        cluster_df = mask_unequal(df, cluster_col, "Noise")

        total_cluster = len(cluster_df)

        cluster_count = len(lst_col(cluster_df, cluster_col, unique=True))

        count_dict = build_col_count_dict(df, cluster_col)
        total_noise = 0
        if "Noise" in list(count_dict.keys()):
            total_noise = count_dict["Noise"]

        cluster_count_lst.append(cluster_count)
        total_noise_lst.append(total_noise)
        grid_graph_lst.append((min_silh, min_simi))

        if report_table_path is not None:
            report_df.at[i, silh_col] = min_silh
            report_df.at[i, simi_col] = min_simi
            report_df.at[i, total_complete_col] = len(df)
            report_df.at[i, cluster_count_col] = cluster_count
            report_df.at[i, total_cluster_pre_col] = total_cluster + total_pruned
            report_df.at[i, total_cluster_post_col] = total_cluster
            report_df.at[i, total_noise_pre_col] = total_noise - total_pruned
            report_df.at[i, total_noise_post_col] = total_noise
            report_df.at[i, total_pruned_nn_col] = nn_pruned
            if constr_matrix is not None:
                report_df.at[i, total_pruned_constr_col] = constr_pruned

            i += 1

    if cluster_lim is not None:
        max_cluster = max(cluster_count_lst)
        for cluster_count in sort_lst(lst_unique(cluster_count_lst)):
            if cluster_count <= cluster_lim:
                max_cluster = cluster_count
            else:
                break

        max_cluster_index_lst = get_lst_val_indices(cluster_count_lst, max_cluster)
        min_noise = min([total_noise_lst[x] for x in max_cluster_index_lst])
        min_noise_index_lst = get_lst_val_indices(total_noise_lst, min_noise)

        best_run = max([x for x in min_noise_index_lst if x in max_cluster_index_lst])
    else:
        min_noise = min(total_noise_lst)
        min_noise_index_lst = get_lst_val_indices(total_noise_lst, min_noise)
        min_cluster = min([cluster_count_lst[x] for x in min_noise_index_lst])
        min_cluster_index_lst = get_lst_val_indices(cluster_count_lst, min_cluster)

        best_run = max([x for x in min_cluster_index_lst if x in min_noise_index_lst])

    if report_table_path is not None:
        report_df.at[best_run, select_col] = True
        report_df[select_col] = report_df[select_col].fillna(False)

    best_min_silh = grid_graph_lst[best_run][0]
    best_min_simi = grid_graph_lst[best_run][1]

    print(f"Best Silhouette={best_min_silh}; Best Similarity={best_min_simi}")

    df = run_grid_graph(
        df,
        result_lst,
        matrix,
        min_min_samples,
        best_min_silh,
        best_min_simi,
        max_dist,
    )

    pre_count = len(lst_col(df, cluster_col, unique=True))

    if constr_matrix is not None:
        if merge_constr_dist is not None:
            if max_constr_dist is None:
                max_constr_dist = build_cutoff_dict(
                    df, constr_matrix, method=constr_method
                )
            df = merge_clusters(
                df, constr_matrix, merge_constr_dist, merge_method=constr_method
            )

    post_merge_count = len(lst_col(df, cluster_col, unique=True))

    result = prune_cluster_members(
        df,
        min_min_samples,
        matrix,
        max_nn_dist=max_nn_dist,
        constr_matrix=constr_matrix,
        max_constr_dist=max_constr_dist,
        constr_method=constr_method,
        min_pdb=min_pdb,
        min_cf=min_cf,
    )

    df = result[0]

    post_prune_count = len(lst_col(df, cluster_col, unique=True))

    print(
        f"Total Clusters={pre_count-1}; Post-Merging={post_merge_count-1}; Post-Pruning={post_prune_count-1}"
    )

    save_table(cluster_table_path, df)

    if sum_table_path is not None:
        sum_df = build_sum_table(df)
        save_table(sum_table_path, sum_df)

    if report_table_path is not None:
        save_table(report_table_path, report_df)

    print("Clustered matrix!")