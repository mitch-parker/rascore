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

from ..scripts.prep_edia import prep_edia
from ..scripts.build_dih_table import build_dih_table
from ..scripts.build_dih_matrix import build_dih_matrix
from ..scripts.build_rmsd_matrix import build_rmsd_matrix
from ..scripts.mask_dih_data import mask_dih_data
from ..scripts.cluster_matrix import cluster_matrix
from ..scripts.classify_matrix import classify_matrix
from ..scripts.build_dist_table import build_dist_table
from ..scripts.write_pymol_script import write_pymol_script


from ..constants.conf import (
    conf_name_dict,
    sw1_name,
    sw2_name,
    nf_name,
    gdp_name,
    gtp_name,
    noise_name,
    loop_resid_dict,
    sw1_gtp_name,
    sw1_gtp_wat_name,
    sw1_gtp_dir_name,
    sw1_gtp_no_name,
    sw1_gtp_dict,
    conf_color_dict,
)
from ..constants.nuc import nuc_class_lst, gtp_atomids
from ..constants.pml import sup_resids, show_resids
from ..functions.cluster import dist_to_dih
from ..functions.col import (
    cluster_col,
    loop_col,
    common_col,
    rama_col,
    rotamer_col,
    pdb_id_col,
    bio_lig_col,
    hb_status_col,
    sw1_col,
    sw2_col,
)
from ..functions.table import lst_col, mask_equal, mask_unequal, make_dict
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
    cluster_report_table_file,
    classify_report_table_file,
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

from ..functions.col import pdb_code_col, nuc_class_col, cluster_col


def cluster_rascore(build_path, out_path=None, name_table_path=None, num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{cluster_str}"

    if name_table_path is None:
        name_dict = conf_name_dict
    else:
        name_df = load_table(name_table_path)

        name_dict = {
            sw1_name: {nf_name: dict(), gdp_name: dict(), gtp_name: dict()},
            sw2_name: {nf_name: dict(), gdp_name: dict(), gtp_name: dict()},
        }
        for index in list(name_df.index.values):
            name_dict[name_df.at[index, "loop"]][name_df.at[index, "nucleotide"]][
                name_df.at[index, "rama"]
            ] = name_df.at[index, "cluster"]

    entry_table_path = get_file_path(entry_table_file, dir_path=build_path)

    sifts_json_path = get_file_path(sifts_json_file, dir_path=build_path)
    edia_json_path = get_file_path(edia_json_file, dir_path=build_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=build_path)

    df = load_table(entry_table_path)

    try:
        pdb_code_lst = lst_col(df, pdb_code_col, unique=True)
        sifts_dict = load_json(sifts_json_path)
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

    for col in [sw1_col, sw2_col]:
        if col in list(df.columns):
            del df[col]

    for loop_name, loop_resids in loop_resid_dict.items():

        loop_result_df = pd.DataFrame()
        loop_sum_df = pd.DataFrame()
        loop_cluster_report_df = pd.DataFrame()
        loop_classify_report_df = pd.DataFrame()

        loop_path = f"{out_path}/{loop_name}"

        for nuc_class in nuc_class_lst:

            loop_nuc_name = f"{loop_name}_{nuc_class}"

            dih_table_path = get_file_path(
                dih_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            fit_table_path = get_file_path(
                fit_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            pred_table_path = get_file_path(
                pred_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            cluster_table_path = get_file_path(
                cluster_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            result_table_path = get_file_path(
                result_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            cluster_report_table_path = get_file_path(
                cluster_report_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )
            classify_report_table_path = get_file_path(
                classify_report_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            sum_table_path = get_file_path(
                sum_table_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            dih_matrix_path = get_file_path(
                dih_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            rmsd_matrix_path = get_file_path(
                rmsd_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            rmsd_json_path = get_file_path(
                rmsd_json_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            dih_fit_matrix_path = get_file_path(
                dih_fit_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            dih_pred_matrix_path = get_file_path(
                dih_pred_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            rmsd_fit_matrix_path = get_file_path(
                rmsd_fit_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            rmsd_pred_matrix_path = get_file_path(
                rmsd_pred_matrix_file, dir_str=loop_nuc_name, dir_path=loop_path
            )

            nuc_df = mask_equal(df, nuc_class_col, nuc_class)

            if len(nuc_df) > 0:

                chi1_resids = None
                if loop_name == sw2_name:
                    chi1_resids = 71

                build_dih_table(
                    df=nuc_df,
                    dih_dict=dih_dict,
                    dih_table_path=dih_table_path,
                    bb_resids=loop_resids,
                    chi1_resids=chi1_resids,
                )
                dih_df = load_table(dih_table_path)

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
                    report_table_path=cluster_report_table_path,
                    max_nn_dist=0.45,
                    constr_matrix=rmsd_fit_matrix,
                    max_constr_dist=1.2,
                    merge_constr_dist=1.2,
                    min_samples_range="3-15",
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
                    report_table_path=classify_report_table_path,
                    sum_table_path=sum_table_path,
                    fit_constr_matrix=rmsd_fit_matrix,
                    pred_constr_matrix=rmsd_pred_matrix,
                    max_nn_dist=0.45,
                    max_constr_dist=1.2,
                )

                result_df = load_table(result_table_path)
                sum_df = load_table(sum_table_path)
                cluster_report_df = load_table(cluster_report_table_path)
                classify_report_df = load_table(classify_report_table_path)

                rename_dict = dict()

                i = 1
                for index in list(sum_df.index.values):
                    cluster = sum_df.at[index, cluster_col]
                    name = cluster
                    if cluster != noise_name:
                        rama = sum_df.at[index, f"{common_col}_{rama_col}"]

                        if loop_name == sw2_name:
                            rama += sum_df.at[index, f"{common_col}_{rotamer_col}"]

                        if rama in list(name_dict[loop_name][nuc_class].keys()):
                            name = name_dict[loop_name][nuc_class][rama]
                        else:
                            name = f"{loop_name}.{nuc_class}-Unknown-{i}"
                            i += 1

                    rename_dict[cluster] = name

                result_df[cluster_col] = result_df[cluster_col].map(rename_dict)
                sum_df[cluster_col] = sum_df[cluster_col].map(rename_dict)
                classify_report_df[cluster_col] = classify_report_df[cluster_col].map(
                    rename_dict
                )

                cluster_df[cluster_col] = cluster_df[pdb_id_col].map(
                    make_dict(
                        lst_col(result_df, pdb_id_col), lst_col(result_df, cluster_col)
                    )
                )

                save_table(cluster_table_path, cluster_df)

                result_df[loop_col] = loop_name
                sum_df[loop_col] = loop_name
                cluster_report_df[loop_col] = loop_name
                classify_report_df[loop_col] = loop_name

                result_df[nuc_class_col] = nuc_class
                sum_df[nuc_class_col] = nuc_class
                cluster_report_df[nuc_class_col] = nuc_class
                classify_report_df[nuc_class_col] = nuc_class

                for index in list(sum_df.index.values):
                    cluster = sum_df.at[index, cluster_col]
                    if cluster == noise_name:
                        sum_df.at[
                            index, cluster_col
                        ] = f"{loop_name}.{nuc_class}-{noise_name}"

                loop_result_df = pd.concat([loop_result_df, result_df], sort=False)
                loop_sum_df = pd.concat([loop_sum_df, sum_df], sort=False)
                loop_cluster_report_df = pd.concat(
                    [loop_cluster_report_df, cluster_report_df], sort=False
                )
                loop_classify_report_df = pd.concat(
                    [loop_classify_report_df, classify_report_df], sort=False
                )

        loop_result_table_path = get_file_path(
            result_table_file, dir_str=loop_name, dir_path=out_path
        )
        loop_sum_table_path = get_file_path(
            sum_table_file, dir_str=loop_name, dir_path=out_path
        )
        loop_cluster_report_table_path = get_file_path(
            cluster_report_table_file, dir_str=loop_name, dir_path=out_path
        )
        loop_classify_report_table_path = get_file_path(
            classify_report_table_file, dir_str=loop_name, dir_path=out_path
        )

        save_table(loop_result_table_path, loop_result_df, fillna="-")
        save_table(loop_sum_table_path, loop_sum_df, fillna="-")
        save_table(loop_cluster_report_table_path, loop_cluster_report_df, fillna="-")
        save_table(loop_classify_report_table_path, loop_classify_report_df, fillna="-")

        cluster_dict = make_dict(
            lst_col(loop_result_df, pdb_id_col), lst_col(loop_result_df, cluster_col)
        )

        df[loop_name] = df[pdb_id_col].map(cluster_dict)

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

        df[hb_status_col] = (
            df[pdb_id_col]
            .map(
                make_dict(lst_col(dist_df, pdb_id_col), lst_col(dist_df, hb_status_col))
            )
            .fillna("")
        )

        df[sw1_name].replace(
            {
                sw1_gtp_name: "",
                sw1_gtp_wat_name: "",
                sw1_gtp_dir_name: "",
                sw1_gtp_no_name: "",
            },
            inplace=True,
        )
        df[sw1_name] += df[hb_status_col].map(str)

        del df[hb_status_col]

    save_table(entry_table_path, df)

    for loop_name, loop_resids in loop_resid_dict.items():

        pymol_pml_path = get_file_path(
            f"{loop_name}_{pymol_pml_file}", dir_path=out_path
        )

        if loop_name == sw1_name:
            stick_resids = [32]
        elif loop_name == sw2_name:
            stick_resids = [71]

        write_pymol_script(
            mask_unequal(df, loop_name, noise_name),
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
            show_resids=show_resids,
        )

    print("Rascore clustering complete!")