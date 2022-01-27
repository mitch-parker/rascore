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
from tqdm import tqdm

from ..functions.cluster import (
    build_cutoff_dict,
    calc_dist_stat,
    order_clusters,
    build_sum_table,
)
from ..functions.table import lst_col, get_val_index_lst, build_col_count_dict
from ..functions.col import (
    cluster_col,
    total_classified_nn_col,
    total_classified_constr_col,
    nn_dist_col,
    constr_dist_col,
    nn_cutoff_col,
    constr_cutoff_col,
    total_col,
)
from ..functions.path import save_table


def classify_matrix(
    cluster_df,
    pred_df,
    fit_matrix,
    pred_matrix,
    result_table_path,
    sum_table_path=None,
    report_table_path=None,
    max_nn_dist=None,
    fit_constr_matrix=None,
    pred_constr_matrix=None,
    max_constr_dist=None,
    constr_method="mean",
    only_save_pred=False,
    reorder_class=True,
    unclass_name="Noise",
):

    if max_nn_dist is None:
        nn_cutoff_dict = build_cutoff_dict(cluster_df, fit_matrix, method="min")
    else:
        nn_cutoff_dict = dict()

    if fit_constr_matrix is not None and pred_constr_matrix is not None:
        if max_constr_dist is None:
            constr_cutoff_dict = build_cutoff_dict(
                cluster_df, fit_constr_matrix, method=constr_method
            )
        else:
            constr_cutoff_dict = dict()

    cluster_lst = lst_col(cluster_df, cluster_col, unique=True, return_str=True)

    report_dict = dict()
    for cluster in cluster_lst:
        report_dict[cluster] = {
            total_classified_nn_col: 0,
            total_classified_constr_col: 0,
        }

    for index in tqdm(
        list(pred_df.index.values),
        desc="Classifying matrix",
        position=0,
        leave=True,
    ):

        pred_lst = list()

        nn_dict = dict()
        constr_dict = dict()

        classified_dict = dict()

        for cluster in cluster_lst:

            add_cluster = True

            index_lst = get_val_index_lst(cluster_df, cluster_col, cluster)

            if len(index_lst) == 1:
                nn_dist = 0
            else:
                if len(pred_df) < 2:
                    dist_lst = pred_matrix[index_lst]
                else:
                    dist_lst = pred_matrix[index, index_lst]
                nn_dist = calc_dist_stat(
                    dist_lst,
                    method="min",
                )

            nn_dict[cluster] = nn_dist

            nn_cutoff = max_nn_dist
            if max_nn_dist is None:
                nn_cutoff = nn_cutoff_dict[cluster]
            else:
                nn_cutoff_dict[cluster] = nn_cutoff

            if nn_dist > nn_cutoff:
                add_cluster = False
            else:
                classified_dict[cluster] = total_classified_nn_col

            if (
                fit_constr_matrix is not None and pred_constr_matrix is not None
            ) and add_cluster is False:
                if len(index_lst) == 1:
                    constr_dist = 0
                else:
                    constr_dist = calc_dist_stat(
                        pred_constr_matrix[index, index_lst],
                        method=constr_method,
                    )
                constr_dict[cluster] = constr_dist

                constr_cutoff = max_constr_dist
                if max_constr_dist is None:
                    constr_cutoff = constr_cutoff_dict[cluster]
                else:
                    constr_cutoff_dict[cluster] = constr_cutoff

                if constr_dist > constr_cutoff:
                    add_cluster = False
                else:
                    classified_dict[cluster] = total_classified_constr_col

            if add_cluster and cluster != unclass_name:
                pred_lst.append(cluster)

        pred_cluster = unclass_name

        if len(pred_lst) == 1:
            pred_cluster = pred_lst[0]
            report_dict[pred_cluster][classified_dict[pred_cluster]] += 1

        pred_df.at[index, cluster_col] = pred_cluster

        if pred_cluster in list(nn_dict.keys()):
            nn_dist = nn_dict[pred_cluster]
        else:
            nn_dist = 0

        pred_df.at[index, nn_dist_col] = nn_dist

        if fit_constr_matrix is not None and pred_constr_matrix is not None:
            if pred_cluster in list(constr_dict.keys()):
                constr_dist = constr_dict[pred_cluster]
            else:
                constr_dist = 0
            pred_df.at[index, constr_dist_col] = constr_dist

    if not only_save_pred:
        pred_df = pd.concat([cluster_df, pred_df], sort=False)
        pred_df = pred_df.reset_index(drop=True)
        if reorder_class:
            pred_df = order_clusters(pred_df)

    save_table(result_table_path, pred_df)

    if sum_table_path is not None:
        sum_df = build_sum_table(pred_df)
        save_table(sum_table_path, sum_df)

    if report_table_path is not None:
        cluster_count_dict = build_col_count_dict(pred_df, cluster_col)
        report_df = pd.DataFrame()
        for i, cluster in enumerate(
            [x for x in cluster_lst if x in list(cluster_count_dict.keys())]
        ):
            report_df.at[i, cluster_col] = cluster
            report_df.at[i, total_col] = cluster_count_dict[cluster]
            report_df.at[i, nn_cutoff_col] = nn_cutoff_dict[cluster]
            report_df.at[i, total_classified_nn_col] = report_dict[cluster][
                total_classified_nn_col
            ]
            if fit_constr_matrix is not None and pred_constr_matrix is not None:
                report_df.at[i, constr_cutoff_col] = constr_cutoff_dict[cluster]
                report_df.at[i, total_classified_constr_col] = report_dict[cluster][
                    total_classified_constr_col
                ]
        save_table(report_table_path, report_df)

    print("Classified matrix!")
