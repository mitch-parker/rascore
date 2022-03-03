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

            if max_nn_dist is None:
                nn_cutoff = nn_cutoff_dict[cluster]
            else:
                nn_cutoff = max_nn_dist
                nn_cutoff_dict[cluster] = nn_cutoff

            if fit_constr_matrix is not None and pred_constr_matrix is not None:
                if max_constr_dist is None:
                    constr_cutoff = constr_cutoff_dict[cluster]
                else:
                    constr_cutoff = max_constr_dist
                    constr_cutoff_dict[cluster] = constr_cutoff

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
            if cluster in list(nn_cutoff_dict.keys()):
                report_df.at[i, nn_cutoff_col] = nn_cutoff_dict[cluster]
            report_df.at[i, total_classified_nn_col] = report_dict[cluster][
                total_classified_nn_col
            ]
            if fit_constr_matrix is not None and pred_constr_matrix is not None:
                if cluster in list(constr_cutoff_dict.keys()):
                    report_df.at[i, constr_cutoff_col] = constr_cutoff_dict[cluster]
                report_df.at[i, total_classified_constr_col] = report_dict[cluster][
                    total_classified_constr_col
                ]
        save_table(report_table_path, report_df)

    print("Classified matrix!")
