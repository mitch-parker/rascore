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

import os
import pandas as pd
from tqdm import tqdm

from ..scripts.annot_lig import annot_lig
from ..scripts.prep_dih import prep_dih
from ..scripts.build_dih_matrix import build_dih_matrix
from ..scripts.build_dih_table import build_dih_table
from ..scripts.build_dist_table import build_dist_table
from ..scripts.classify_matrix import classify_matrix
from ..scripts.write_pymol_script import write_pymol_script

from ..functions.table import (
    mask_equal,
    lst_to_str,
    lst_col,
    type_lst,
    merge_tables,
    make_dict,
    fix_col,
    core_path_col,
    chainid_col,
    modelid_col,
    pdb_id_col,
)
from ..functions.path import (
    load_matrix,
    save_table,
    get_file_name,
    get_file_path,
    load_lst,
    load_table,
    get_neighbor_path,
    rascore_str,
    classify_str,
    cluster_str,
    pipelines_str,
    data_str,
)
from ..functions.cluster import build_sum_table
from ..functions.coord import load_coord, get_modelid, get_chainid
from ..functions.lig import lig_col_lst
from ..functions.col import (
    id_col,
    nuc_class_col,
    bio_lig_col,
    cluster_col,
    hb_status_col,
    complete_col,
)
from ..functions.file import (
    cluster_table_file,
    dih_fit_matrix_file,
    result_table_file,
    sum_table_file,
    pymol_pml_file,
    pred_matrix_file,
)
from ..constants.nuc import (
    gtp_name,
    nuc_class_dict,
    nuc_class_dict,
    gtp_atomids,
    nuc_class_lst,
)
from ..constants.pharm import (
    pharm_site_dict,    
)
from ..constants.conf import (
    y32_name,
    y71_name,
    sw1_name,
    sw2_name,
    in_name,
    out_name,
    noise_name,
    outlier_name,
    disorder_name,
    conf_color_dict
)
from ..constants.pml import sup_resids, show_resids


def classify_rascore(file_paths, out_path=None, dih_dict=None, 
                    sw1_resids="25-40", sw2_resids="56-76", 
                    y32_resid=32, y71_resid=71, 
                    g12_resid=12, v9_resid=9,
                    y32_dist=10.5, y71_dist=8.75,
                    num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{classify_str}"

    cluster_path = f"{get_neighbor_path(__file__,pipelines_str,data_str)}/{rascore_str}_{cluster_str}"

    file_path_lst = type_lst(file_paths)

    df = pd.DataFrame()
    for file_path in tqdm(
        file_path_lst,
        desc="Loading Files",
        position=0,
        leave=True,
    ):
        if type(file_path) == pd.DataFrame:
            temp_df = file_path.copy(deep=True)
            temp_df_col_lst = list(temp_df.columns)
            if modelid_col not in temp_df_col_lst:
                temp_df[modelid_col] = 0
        else:
            temp_df = pd.DataFrame()
            load_path_lst = list()

            if ".pdb" not in file_path and ".cif" not in file_path:
                temp_df = load_table(file_path)
                temp_df_col_lst = list(temp_df.columns)
                if core_path_col in temp_df_col_lst and chainid_col in temp_df_col_lst:
                    if modelid_col not in temp_df_col_lst:
                        temp_df[modelid_col] = 0
                else:
                    load_path_lst = load_lst(file_path)
                    temp_df = pd.DataFrame()
            else:
                load_path_lst = list([file_path])

            if len(temp_df) == 0 and len(load_path_lst) > 0:
                i = 0
                for load_path in tqdm(
                    load_path_lst,
                    desc="Loading Files",
                    position=1,
                    leave=True,
                ):
                    structure = load_coord(load_path)
                    for model in structure:
                        modelid = get_modelid(model)
                        for chain in model:
                            chainid = get_chainid(chain)
                            temp_df.at[i, core_path_col] = load_path
                            temp_df.at[i, modelid_col] = modelid
                            temp_df.at[i, chainid_col] = chainid
                            i += 1

        df = pd.concat([df, temp_df], sort=False)

    if len(df) == 0:
        print('No structures to classify.')
    else:
        df = df.reset_index(drop=True)

        df_col_lst = list(df.columns)

        for col in [y32_name, y71_name, sw1_name, sw2_name]:
            if col in df_col_lst:
                del df[col]

        missing_col_lst = [
            x for x in [core_path_col, modelid_col, chainid_col] if x not in df_col_lst
        ]
        if len(missing_col_lst) > 0:
            print(f"Input table missing columns: {lst_to_str(missing_col_lst)}")

        else:
            coord_path_lst = lst_col(df, core_path_col)

            if id_col not in df_col_lst and pdb_id_col not in df_col_lst:
                for index in list(df.index.values):
                    df.at[index, id_col] = get_file_name(df.at[index, core_path_col])

            if len([x for x in lig_col_lst if x not in df_col_lst]) > 0:
                df = annot_lig(
                    df=df,
                    site_dict=pharm_site_dict,
                    num_cpu=num_cpu,
                )

            if nuc_class_col not in df_col_lst:
                df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

            if dih_dict is None:
                dih_dict = prep_dih(coord_path_lst, num_cpu=num_cpu)

            dist_df = build_dist_table(
                df,
                x_resids=[y32_resid, y71_resid],
                y_resids=[g12_resid, v9_resid],
                x_atomids=['OH','OH'],
                y_atomids=['CA','CA'],
                atom_dist_col_lst=[y32_name, y71_name]
            )

            dist_dict = {y32_name: y32_dist, y71_name: y71_dist}

            for resid_name in [y32_name, y71_name]:
                dist_df[resid_name] = dist_df[resid_name].map(str)

            for index in list(dist_df.index.values):

                nuc_class = dist_df.at[index, nuc_class_col]

                for resid_name in [y32_name, y71_name]:

                    atom_dist = float(dist_df.at[index, resid_name])

                    if atom_dist == 999.00:
                        group_name = disorder_name
                    else:
                        group_name = resid_name
                        if atom_dist <= dist_dict[resid_name]:
                            group_name += in_name
                        else:
                            group_name += out_name

                    dist_df.at[index, resid_name] = group_name

            for resid_name in [y32_name, y71_name]:

                merge_df = dist_df.loc[
                    :, [core_path_col, modelid_col, chainid_col, resid_name]
                ]

                df = merge_tables(df, merge_df)

            for resid_name in [y32_name, y71_name]:

                if resid_name == y32_name:
                    loop_name = sw1_name
                    loop_resids = sw1_resids
                elif resid_name == y71_name:
                    loop_name = sw2_name
                    loop_resids = sw2_resids

                cluster_loop_path = f"{cluster_path}/{loop_name}"
                classify_loop_path = f"{out_path}/{loop_name}"

                loop_result_df = pd.DataFrame()
                loop_sum_df = pd.DataFrame()

                for group_name in lst_col(df, resid_name, unique=True):

                    group_df = mask_equal(df, resid_name, group_name)

                    for nuc_class in nuc_class_lst:

                        group_nuc_df = mask_equal(group_df, nuc_class_col, nuc_class)

                        group_nuc_name = f"{group_name}.{nuc_class}"

                        print(
                                f"Classifying {loop_name} conformations in {group_name}.{nuc_class} structures."
                            )

                        cluster_table_path = get_file_path(
                            cluster_table_file,
                            dir_str=group_nuc_name,
                            dir_path=cluster_loop_path,
                        )

                        fit_matrix_path = get_file_path(
                            dih_fit_matrix_file,
                            dir_str=group_nuc_name,
                            dir_path=cluster_loop_path,
                        )

                        result_table_path = get_file_path(
                            result_table_file,
                            dir_str=group_nuc_name,
                            dir_path=classify_loop_path,
                        )

                        pred_matrix_path = get_file_path(
                            pred_matrix_file,
                            dir_str=group_nuc_name,
                            dir_path=classify_loop_path,
                        )

                        if len(group_nuc_df) > 0:

                            dih_df = build_dih_table(
                                    df=group_nuc_df,
                                    dih_dict=dih_dict,
                                    bb_resids=loop_resids,
                                )

                            if disorder_name not in group_nuc_name:

                                cluster_df = load_table(cluster_table_path)

                                build_dih_matrix(
                                    fit_df=cluster_df,
                                    pred_df=dih_df,
                                    max_norm_path=pred_matrix_path,
                                )

                                fit_matrix = load_matrix(fit_matrix_path)
                                pred_matrix = load_matrix(pred_matrix_path)

                                classify_matrix(
                                    cluster_df=cluster_df,
                                    pred_df=dih_df,
                                    fit_matrix=fit_matrix,
                                    pred_matrix=pred_matrix,
                                    result_table_path=result_table_path,
                                    max_nn_dist=0.45,
                                    only_save_pred=True,
                                    reorder_class=False,
                                )
                            else:
                                result_df = dih_df.copy(deep=True)
                                result_df[cluster_col] = noise_name
                                save_table(result_table_path, result_df)

                            result_df = load_table(result_table_path)

                            loop_result_df = pd.concat(
                                [loop_result_df, result_df], sort=False
                            )
                       
                loop_result_table_path = get_file_path(
                    result_table_file, dir_str=loop_name, dir_path=out_path
                )
                loop_sum_table_path = get_file_path(
                    sum_table_file, dir_str=loop_name, dir_path=out_path
                )

                if pdb_id_col in list(loop_result_df.columns):

                    cluster_result_df = load_table(
                        get_file_path(
                            result_table_file,
                            dir_str=loop_name,
                            dir_path=cluster_path,
                        )
                    )

                    cluster_dict = make_dict(
                        lst_col(cluster_result_df, pdb_id_col),
                        lst_col(cluster_result_df, cluster_col),
                    )

                    result_dict = make_dict(
                        lst_col(loop_result_df, pdb_id_col),
                        lst_col(loop_result_df, cluster_col),
                    )

                    for pdb_id, cluster in cluster_dict.items():
                        if pdb_id in list(result_dict.keys()):
                            result_dict[pdb_id] = cluster

                    loop_result_df[cluster_col] = loop_result_df[pdb_id_col].map(
                        result_dict
                    )

                loop_result_df = loop_result_df.reset_index(drop=True)
            
                for index in list(loop_result_df.index.values):
                    cluster = loop_result_df.at[index, cluster_col]
                    if noise_name in cluster:
                        complete =  str(loop_result_df.at[index, complete_col])
                        if complete == str(True):
                            cluster = outlier_name
                        elif complete == str(False):
                            cluster = disorder_name
                        loop_result_df.at[index,cluster_col] = cluster

                loop_sum_df = build_sum_table(loop_result_df)

                save_table(loop_result_table_path, loop_result_df)
                save_table(loop_sum_table_path, loop_sum_df)

                loop_result_df = loop_result_df.rename(columns={cluster_col: loop_name})

                loop_result_df = loop_result_df.loc[:,[core_path_col, modelid_col, chainid_col, loop_name]]

                df = merge_tables(df, loop_result_df)

            if gtp_name in lst_col(df, nuc_class_col, unique=True):
                dist_df = build_dist_table(
                    mask_equal(df, nuc_class_col, gtp_name),
                    x_resids=[y32_resid],
                    y_resids=[bio_lig_col],
                    x_atomids=["OH"],
                    y_atomids=[gtp_atomids],
                    hb_status_col_lst=[hb_status_col],
                    check_hb=True,
                )

                merge_df = dist_df.loc[
                    :, [core_path_col, modelid_col, chainid_col, hb_status_col]
                ]

                df = merge_tables(df, merge_df)

            result_table_path = get_file_path(result_table_file, dir_path=out_path)

            for col in list(df.columns):
                df = fix_col(df, col)

            save_table(result_table_path, df)

            for loop_name in [sw1_name,sw2_name]:

                pymol_pml_path = get_file_path(
                    f"{loop_name}_{pymol_pml_file}", dir_path=out_path
                )

                if loop_name == sw1_name:
                    loop_resids = sw1_resids
                    stick_resids = [y32_resid]
                elif loop_name == sw2_name:
                    loop_resids = sw2_resids
                    stick_resids = [y71_resid]

                write_pymol_script(
                    df,
                    pymol_pml_path,
                    group_col=loop_name,
                    stick_resids=stick_resids,
                    loop_resids=loop_resids,
                    style_ribbon=True,
                    thick_bb=False,
                    show_bio=True,
                    sup_group=True,
                    color_palette=conf_color_dict[loop_name],
                    sup_resids=sup_resids,
                    show_resids=show_resids,
                )

    print("Rascore classification complete!")