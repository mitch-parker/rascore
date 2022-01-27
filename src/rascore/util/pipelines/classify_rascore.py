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
from ..functions.coord import load_coord, get_modelid, get_chainid
from ..functions.lig import lig_col_lst
from ..functions.col import (
    id_col,
    pharm_class_col,
    pharm_lig_site_col,
    nuc_class_col,
    bio_lig_col,
    loop_col,
    cluster_col,
    hb_status_col,
)
from ..functions.file import (
    cluster_table_file,
    dih_fit_matrix_file,
    result_table_file,
    sum_table_file,
    classify_report_table_file,
    pymol_pml_file,
    pred_matrix_file,
)
from ..constants.nuc import (
    gtp_name,
    nuc_class_dict,
    nuc_class_lst,
    nuc_class_dict,
    gtp_atomids,
)
from ..constants.pharm import (
    pharm_site_dict,
    pharm_match_dict,
    sp2_name,
    sp12_name,
    none_pharm_name,
    other_pharm_name,
)
from ..constants.conf import (
    loop_resid_dict,
    sw1_name,
    sw2_name,
    sw1_gtp_name,
    sw1_gtp_dir_name,
    sw1_gtp_wat_name,
    sw1_gtp_no_name,
    sw1_gtp_dict,
    conf_color_dict,
)
from ..constants.pml import sup_resids


def classify_rascore(coord_paths, out_path=None, dih_dict=None, num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{classify_str}"

    cluster_path = f"{get_neighbor_path(__file__,pipelines_str,data_str)}/{rascore_str}_{cluster_str}"

    coord_path_lst = type_lst(coord_paths)

    if type(coord_path_lst[0]) == pd.DataFrame:
        df = coord_path_lst[0]
    else:
        df = pd.DataFrame()
        if ".pdb" not in coord_path_lst[0] and ".cif" not in coord_path_lst[0]:
            df = load_table(coord_path_lst[0])
            df_col_lst = list(df.columns)
            if core_path_col in df_col_lst and chainid_col in df_col_lst:
                if modelid_col not in df_col_lst:
                    df[modelid_col] = 0
            else:
                if core_path_col in df_col_lst and chainid_col not in df_col_lst:
                    coord_path_lst = lst_col(df, core_path_col)
                else:
                    coord_path_lst = load_lst(coord_path_lst[0])
                df = pd.DataFrame()

        if len(df) == 0:

            i = 0
            for coord_path in tqdm(
                coord_path_lst,
                desc="Loading coordinates",
                position=0,
                leave=True,
            ):
                structure = load_coord(coord_path)
                for model in structure:
                    modelid = get_modelid(model)
                    for chain in model:
                        chainid = get_chainid(chain)
                        df.at[i, core_path_col] = coord_path
                        df.at[i, modelid_col] = modelid
                        df.at[i, chainid_col] = chainid
                        i += 1

    if len(df) > 0:
        df_col_lst = list(df.columns)

        if id_col not in df_col_lst and pdb_id_col not in df_col_lst:
            for index in list(df.index.values):
                df.at[index, id_col] = get_file_name(df.at[index, core_path_col])

        missing_col_lst = [
            x for x in [core_path_col, modelid_col, chainid_col] if x not in df_col_lst
        ]
        if len(missing_col_lst) > 0:
            print(f"Input table missing columns: {lst_to_str(missing_col_lst)}")

        else:
            if len([x for x in lig_col_lst if x not in df_col_lst]) > 0:
                df = annot_lig(
                    df=df,
                    site_dict=pharm_site_dict,
                    match_dict=pharm_match_dict,
                    num_cpu=num_cpu,
                )

            if pharm_class_col not in df_col_lst:
                for index in list(df.index.values):
                    pharm_class = df.at[index, pharm_lig_site_col]
                    if pharm_class not in [sp2_name, sp12_name, none_pharm_name]:
                        pharm_class = other_pharm_name
                    df.at[index, pharm_class_col] = pharm_class

            if nuc_class_col not in df_col_lst:
                df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

            if dih_dict is None:
                dih_dict = prep_dih(coord_paths, num_cpu=num_cpu)

            for loop_name, loop_resids in loop_resid_dict.items():

                cluster_loop_path = f"{cluster_path}/{loop_name}"
                classify_loop_path = f"{out_path}/{loop_name}"

                loop_result_df = pd.DataFrame()
                loop_sum_df = pd.DataFrame()
                loop_classify_report_df = pd.DataFrame()

                for nuc_class in nuc_class_lst:

                    print(
                        f"Classifying {loop_name} conformations in {nuc_class} structures."
                    )

                    loop_nuc_name = f"{loop_name}_{nuc_class}"

                    cluster_table_path = get_file_path(
                        cluster_table_file,
                        dir_str=loop_nuc_name,
                        dir_path=cluster_loop_path,
                    )

                    fit_matrix_path = get_file_path(
                        dih_fit_matrix_file,
                        dir_str=loop_nuc_name,
                        dir_path=cluster_loop_path,
                    )

                    result_table_path = get_file_path(
                        result_table_file,
                        dir_str=loop_nuc_name,
                        dir_path=classify_loop_path,
                    )

                    sum_table_path = get_file_path(
                        sum_table_file,
                        dir_str=loop_nuc_name,
                        dir_path=classify_loop_path,
                    )

                    classify_report_table_path = get_file_path(
                        classify_report_table_file,
                        dir_str=loop_nuc_name,
                        dir_path=classify_loop_path,
                    )

                    pred_matrix_path = get_file_path(
                        pred_matrix_file,
                        dir_str=loop_nuc_name,
                        dir_path=classify_loop_path,
                    )

                    nuc_df = mask_equal(df, nuc_class_col, nuc_class)

                    if len(nuc_df) > 0:

                        chi1_resids = None
                        if loop_name == sw2_name:
                            chi1_resids = 71

                        dih_df = build_dih_table(
                            df=nuc_df,
                            dih_dict=dih_dict,
                            bb_resids=loop_resids,
                            chi1_resids=chi1_resids,
                        )

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
                            sum_table_path=sum_table_path,
                            report_table_path=classify_report_table_path,
                            max_nn_dist=0.45,
                            only_save_pred=True,
                            reorder_class=False,
                        )

                        result_df = load_table(result_table_path)
                        sum_df = load_table(sum_table_path)
                        classify_report_df = load_table(classify_report_table_path)

                        result_df[loop_col] = loop_name
                        sum_df[loop_col] = loop_name
                        classify_report_df[loop_col] = loop_name

                        result_df[nuc_class_col] = nuc_class
                        sum_df[nuc_class_col] = nuc_class
                        classify_report_df[nuc_class_col] = nuc_class

                        loop_result_df = pd.concat(
                            [loop_result_df, result_df], sort=False
                        )
                        loop_sum_df = pd.concat([loop_sum_df, sum_df], sort=False)
                        loop_classify_report_df = pd.concat(
                            [loop_classify_report_df, classify_report_df], sort=False
                        )

                loop_result_table_path = get_file_path(
                    result_table_file, dir_str=loop_name, dir_path=out_path
                )
                loop_sum_table_path = get_file_path(
                    sum_table_file, dir_str=loop_name, dir_path=out_path
                )
                loop_classify_report_table_path = get_file_path(
                    classify_report_table_file,
                    dir_str=loop_name,
                    dir_path=out_path,
                )

                save_table(loop_result_table_path, loop_result_df)
                save_table(loop_sum_table_path, loop_sum_df)
                save_table(loop_classify_report_table_path, loop_classify_report_df)

                loop_result_df = loop_result_df.rename(columns={cluster_col: loop_name})

                loop_result_df = loop_result_df.loc[
                    :, [core_path_col, modelid_col, chainid_col, loop_name]
                ]

                df = merge_tables(df, loop_result_df)

            if sw1_gtp_name in lst_col(df, sw1_name, unique=True):
                dist_df = build_dist_table(
                    mask_equal(df, sw1_name, sw1_gtp_name),
                    x_resids=[32],
                    y_resids=[bio_lig_col],
                    x_atomids=["OH"],
                    y_atomids=[gtp_atomids],
                    hb_status_col_lst=[hb_status_col],
                    check_hb=True,
                )

                dist_df[hb_status_col] = dist_df[hb_status_col].map(sw1_gtp_dict)

                dist_df = dist_df.loc[
                    :, [core_path_col, modelid_col, chainid_col, hb_status_col]
                ]

                df = merge_tables(df, dist_df)

                df[sw1_name].replace(
                    {
                        sw1_gtp_name: "",
                        sw1_gtp_wat_name: "",
                        sw1_gtp_dir_name: "",
                        sw1_gtp_no_name: "",
                    },
                    inplace=True,
                )
                df[sw1_name] += df[hb_status_col].fillna("").map(str)

                del df[hb_status_col]

            result_table_path = get_file_path(result_table_file, dir_path=out_path)

            save_table(result_table_path, df)

            for loop_name, loop_resids in loop_resid_dict.items():

                pymol_pml_path = get_file_path(
                    f"{loop_name}_{pymol_pml_file}", dir_path=out_path
                )

                if loop_name == sw1_name:
                    stick_resids = [32]
                elif loop_name == sw2_name:
                    stick_resids = [71]

                write_pymol_script(
                    df,
                    pymol_pml_path,
                    group_col=loop_name,
                    stick_resids=stick_resids,
                    loop_resids=loop_resids,
                    style_ribbon=True,
                    thick_bb=False,
                    show_bio=True,
                    color_palette=conf_color_dict[loop_name],
                    sup_group=True,
                    sup_resids=sup_resids,
                    show_resids="1-166",
                )

    print("Rascore classification complete!")