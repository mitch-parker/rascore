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

from ..scripts.annot_lig import annot_lig
from ..scripts.prep_dih import prep_dih
from ..scripts.prep_edia import prep_edia
from ..scripts.build_dih_table import build_dih_table
from ..scripts.build_dih_matrix import build_dih_matrix
from ..scripts.build_rmsd_matrix import build_rmsd_matrix
from ..scripts.mask_dih_data import mask_dih_data
from ..scripts.cluster_matrix import cluster_matrix
from ..scripts.classify_matrix import classify_matrix
from ..scripts.build_dist_table import build_dist_table
from ..scripts.write_pymol_script import write_pymol_script

from ..functions.cluster import build_sum_table

from ..constants.conf import (
    conf_name_dict,
    y32_name,
    y71_name,
    sw1_name,
    sw2_name,
    in_name,
    out_name,
    outlier_name,
    disorder_name,
    gtp_name,
    noise_name,
    conf_color_dict
)
from ..constants.nuc import nuc_class_dict, nuc_class_lst
from ..constants.pml import sup_resids, show_resids

from ..functions.table import get_col_most_common, lst_col, mask_equal, make_dict
from ..functions.path import (
    load_table,
    save_table,
    load_matrix,
    delete_path,
    load_table,
    load_json,
    get_file_path,
    delete_path,
    rascore_str,
    cluster_str,
)
from ..functions.file import (
    entry_table_file,
    sifts_json_file,
    edia_json_file,
    dih_json_file,
    dih_table_file,
    fit_table_file,
    pred_table_file,
    cluster_table_file,
    result_table_file,
    sum_table_file,
    dih_matrix_file,
    rmsd_matrix_file,
    rmsd_json_file,
    dih_fit_matrix_file,
    dih_pred_matrix_file,
    rmsd_fit_matrix_file,
    rmsd_pred_matrix_file,
    pymol_pml_file,
)

from ..functions.lig import lig_col_lst
from ..functions.col import (
    pdb_code_col, 
    complete_col,
    rama_col,
    pdb_id_col,
    bio_lig_col,
    nuc_class_col, 
    cluster_col, 
    core_path_col,)


def cluster_rascore(build_path, out_path=None, name_dict=None,   
                    sw1_resids="25-40", sw2_resids="56-76", 
                    g12_resid=12, v9_resid=9,
                    y32_resid=32, y71_resid=71, 
                    y32_dist=10.5, y71_dist=8.75,
                    num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{cluster_str}"

    if name_dict is None:
        name_dict = conf_name_dict

    entry_table_path = get_file_path(entry_table_file, dir_path=build_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=build_path)

    df = load_table(entry_table_path)

    df_col_lst = list(df.columns)

    for col in [y32_name, y71_name, sw1_name, sw2_name]:
        if col in df_col_lst:
            del df[col]

    if pdb_code_col in df_col_lst:
        try:
            edia_json_path = get_file_path(edia_json_file, dir_path=build_path)
            sifts_json_path = get_file_path(sifts_json_file, dir_path=build_path)

            sifts_dict = load_json(sifts_json_path)

            pdb_code_lst = lst_col(df, pdb_code_col, unique=True)
            prep_edia(
                pdb_codes=pdb_code_lst,
                edia_dir=build_path,
                sifts_dict=sifts_dict,
                edia_json_path=edia_json_path,
                num_cpu=num_cpu,
            )
        except:
            delete_path(edia_json_path)

    dih_dict = load_json(dih_json_path)
    edia_dict = load_json(edia_json_path)

    if dih_dict is None:
        prep_dih(lst_col(df,core_path_col),dih_json_path=dih_json_path,num_cpu=num_cpu)
        dih_dict = load_json(dih_json_path)

    if len([x for x in lig_col_lst if x not in df_col_lst]) > 0:
        df = annot_lig(
            df=df,
            num_cpu=num_cpu,
        )

    if nuc_class_col not in df_col_lst:
        df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

    dist_df = build_dist_table(
            df,
            x_resids=[y32_resid, y71_resid],
            y_resids=[g12_resid, v9_resid],
            x_atomids=['OH','OH'],
            y_atomids=['CA','CA'],
            atom_dist_col_lst=[y32_name, y71_name]
        )

    dist_dict = {y32_name: y32_dist, y71_name: y71_dist}

    for index in list(dist_df.index.values):

        nuc_class = dist_df.at[index, nuc_class_col]

        for resid_name in [y32_name, y71_name]:

            atom_dist = dist_df.at[index, resid_name]

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

        df[resid_name] = (
            df[pdb_id_col]
            .map(
                make_dict(lst_col(dist_df, pdb_id_col), lst_col(dist_df, resid_name))
            )
        )

    for resid_name in [y32_name, y71_name]:

        if resid_name == y32_name:
            loop_name = sw1_name
            loop_resids = sw1_resids
        elif resid_name == y71_name:
            loop_name = sw2_name
            loop_resids = sw2_resids

        loop_result_df = pd.DataFrame()
        loop_sum_df = pd.DataFrame()

        loop_path = f"{out_path}/{loop_name}"

        total_complete = 0
        total_clustered = 0
        
        for group_name in lst_col(df, resid_name, unique=True):

            group_df = mask_equal(df, resid_name, group_name)

            for nuc_class in nuc_class_lst:

                group_nuc_df = mask_equal(group_df, nuc_class_col, nuc_class)

                group_nuc_name = f"{group_name}.{nuc_class}"

                dih_table_path = get_file_path(
                    dih_table_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                fit_table_path = get_file_path(
                    fit_table_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                pred_table_path = get_file_path(
                    pred_table_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                cluster_table_path = get_file_path(
                    cluster_table_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                result_table_path = get_file_path(
                    result_table_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                dih_matrix_path = get_file_path(
                    dih_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                rmsd_matrix_path = get_file_path(
                    rmsd_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                rmsd_json_path = get_file_path(
                    rmsd_json_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                dih_fit_matrix_path = get_file_path(
                    dih_fit_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                dih_pred_matrix_path = get_file_path(
                    dih_pred_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                rmsd_fit_matrix_path = get_file_path(
                    rmsd_fit_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                rmsd_pred_matrix_path = get_file_path(
                    rmsd_pred_matrix_file, dir_str=group_nuc_name, dir_path=loop_path
                )

                build_dih_table(
                    df=group_nuc_df,
                    dih_dict=dih_dict,
                    dih_table_path=dih_table_path,
                    bb_resids=loop_resids,
                )
                dih_df = load_table(dih_table_path)

                if len(dih_df) > 15 and disorder_name not in group_nuc_name:

                    rmsd_dict = load_json(rmsd_json_path)
                    build_rmsd_matrix(
                        fit_df=dih_df,
                        rmsd_matrix_path=rmsd_matrix_path,
                        sup_resids=sup_resids,
                        rmsd_resids=loop_resids,
                        rmsd_atomids="CA",
                        pair_aln=False,
                        rmsd_dict=rmsd_dict,
                        rmsd_json_path=rmsd_json_path,
                    )

                    build_dih_matrix(
                        fit_df=dih_df,
                        max_norm_path=dih_matrix_path,
                    )

                    dih_matrix = load_matrix(dih_matrix_path)
                    rmsd_matrix = load_matrix(rmsd_matrix_path)

                    mask_dih_data(
                        df=dih_df,
                        matrix=dih_matrix,
                        fit_table_path=fit_table_path,
                        fit_matrix_path=dih_fit_matrix_path,
                        pred_table_path=pred_table_path,
                        pred_matrix_path=dih_pred_matrix_path,
                        edia_dict=edia_dict,
                        edia_min=0.4,
                        edia_atomids="O",
                    )
                    mask_dih_data(
                        df=dih_df,
                        matrix=rmsd_matrix,
                        fit_table_path=fit_table_path,
                        fit_matrix_path=rmsd_fit_matrix_path,
                        pred_table_path=pred_table_path,
                        pred_matrix_path=rmsd_pred_matrix_path,
                        edia_dict=edia_dict,
                        edia_min=0.4,
                        edia_atomids="O",
                    )

                    fit_df = load_table(fit_table_path)
                    dih_fit_matrix = load_matrix(dih_fit_matrix_path)
                    rmsd_fit_matrix = load_matrix(rmsd_fit_matrix_path)

                    cluster_matrix(
                        df=fit_df,
                        matrix=dih_fit_matrix,
                        cluster_table_path=cluster_table_path,
                        max_nn_dist=0.45,
                        constr_matrix=rmsd_fit_matrix,
                        max_constr_dist=1.2,
                        merge_constr_dist=1.2,
                        min_samples_range="7-15",
                        min_min_samples=7,
                        min_pdb=5,
                    )

                    cluster_df = load_table(cluster_table_path)

                    pred_df = load_table(pred_table_path)

                    dih_pred_matrix = load_matrix(dih_pred_matrix_path)
                    rmsd_pred_matrix = load_matrix(rmsd_pred_matrix_path)

                    classify_matrix(
                        cluster_df=cluster_df,
                        pred_df=pred_df,
                        fit_matrix=dih_fit_matrix,
                        pred_matrix=dih_pred_matrix,
                        result_table_path=result_table_path,
                        fit_constr_matrix=rmsd_fit_matrix,
                        pred_constr_matrix=rmsd_pred_matrix,
                        max_nn_dist=0.45,
                        max_constr_dist=1.2,
                    )
                else:
                    cluster_df = dih_df.copy(deep=True)
                    cluster_df[cluster_col] = noise_name

                    result_df = cluster_df.copy(deep=True)

                    save_table(cluster_table_path, cluster_df)
                    save_table(result_table_path, result_df)

                cluster_df = load_table(cluster_table_path)
                result_df = load_table(result_table_path)

                rename_dict = dict()

                i = 1
                for cluster in lst_col(result_df,cluster_col,unique=True):
                    name = cluster
                    if cluster != noise_name:
                        rama = [x for x in get_col_most_common(mask_equal(result_df,cluster_col,cluster), rama_col) if "-" not in x][0]
                        name = f"{group_nuc_name}-Unknown-{i}"
                        i += 1
                        if group_nuc_name in list(name_dict.keys()):
                            if rama in list(name_dict[group_nuc_name].keys()):
                                name = name_dict[group_nuc_name][rama]
                                i -= 1                        

                    rename_dict[cluster] = name

                result_df[cluster_col] = result_df[cluster_col].map(rename_dict)

                cluster_df[cluster_col] = cluster_df[pdb_id_col].map(
                    make_dict(
                        lst_col(result_df, pdb_id_col), lst_col(result_df, cluster_col)
                    )
                )

                save_table(cluster_table_path, cluster_df)

                for index in list(result_df.index.values):
                    cluster = result_df.at[index, cluster_col]
                    if noise_name in cluster:
                        complete =  str(result_df.at[index, complete_col])
                        if complete == str(True):
                            cluster = outlier_name
                        elif complete == str(False):
                            cluster = disorder_name
                        result_df.at[index,cluster_col] = cluster

                loop_result_df = pd.concat([loop_result_df, result_df], sort=False)

                total_complete += len(cluster_df)
                total_clustered += len(mask_equal(cluster_df, cluster_col, [x for x in lst_col(cluster_df, cluster_col) if noise_name not in x]))

        loop_result_table_path = get_file_path(
            result_table_file, dir_str=loop_name, dir_path=out_path
        )
        loop_sum_table_path = get_file_path(
            sum_table_file, dir_str=loop_name, dir_path=out_path
        )

        loop_sum_df = build_sum_table(loop_result_df)

        save_table(loop_result_table_path, loop_result_df, fillna="-")
        save_table(loop_sum_table_path, loop_sum_df, fillna="-")

        cluster_dict = make_dict(
            lst_col(loop_result_df, pdb_id_col), lst_col(loop_result_df, cluster_col)
        )

        df[loop_name] = df[pdb_id_col].map(cluster_dict)

        final_spatial = len(mask_equal(loop_result_df, resid_name, [x for x in lst_col(loop_result_df, resid_name) if disorder_name not in x]))
        final_conformation = len(mask_equal(loop_result_df, cluster_col, [x for x in lst_col(loop_result_df, cluster_col) if outlier_name not in x and disorder_name not in x]))

        print(f"Results for {loop_name}: Total Well Modeled = {total_complete}; Total Clustered = {total_clustered}; Final with Spatial Labels = {final_spatial}; Final with Conformational Labels = {final_conformation}")

    save_table(entry_table_path, df)

    for loop_name in [sw1_name, sw2_name]:

        if loop_name == sw1_name:
            loop_resids = sw1_resids
        elif loop_name == sw2_name:
            loop_resids = sw2_resids

        pymol_pml_path = get_file_path(
            f"{loop_name}_{pymol_pml_file}", dir_path=out_path
        )

        if loop_name == sw1_name:
            stick_resids = [32]
        elif loop_name == sw2_name:
            stick_resids = [71]

        write_pymol_script(
            mask_equal(df, loop_name, [x for x in lst_col(df, loop_name, unique=True) if noise_name not in x]),
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

    print("Rascore clustering complete!")