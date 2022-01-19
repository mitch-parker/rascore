# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from .scripts import *
from .constants import *


def cluster_rascore(data_path, out_path=None, name_table_path=None, num_cpu=1):

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

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)

    dih_json_path = get_file_path(dih_json_file, dir_path=data_path)
    edia_json_path = get_file_path(edia_json_file, dir_path=data_path)

    df = load_table(entry_table_path)
    dih_dict = load_json(dih_json_path)
    edia_dict = load_json(edia_json_path)

    for loop_name, loop_resids in loop_resid_dict.items():

        loop_result_df = pd.DataFrame()
        loop_sum_df = pd.DataFrame()
        loop_cluster_report_df = pd.DataFrame()
        loop_classify_report_df = pd.DataFrame()

        for nuc_class in nuc_class_lst:

            run_name = f"{loop_name}_{nuc_class}"

            dih_table_path = get_file_path(
                dih_table_file, dir_str=run_name, dir_path=data_path
            )

            included_table_path = get_file_path(
                included_table_file, dir_str=run_name, dir_path=data_path
            )

            removed_table_path = get_file_path(
                removed_table_file, dir_str=run_name, dir_path=data_path
            )

            cluster_table_path = get_file_path(
                cluster_table_file, dir_str=run_name, dir_path=data_path
            )

            result_table_path = get_file_path(
                result_table_file, dir_str=run_name, dir_path=data_path
            )

            cluster_report_table_path = get_file_path(
                cluster_report_table_file, dir_str=run_name, dir_path=data_path
            )
            classify_report_table_path = get_file_path(
                classify_report_table_file, dir_str=run_name, dir_path=data_path
            )

            sum_table_path = get_file_path(
                sum_table_file, dir_str=run_name, dir_path=data_path
            )

            dih_matrix_path = get_file_path(
                dih_matrix_file, dir_str=run_name, dir_path=data_path
            )

            rmsd_matrix_path = get_file_path(
                rmsd_matrix_file, dir_str=run_name, dir_path=data_path
            )

            rmsd_json_path = get_file_path(
                rmsd_json_file, dir_str=run_name, dir_path=data_path
            )

            dih_fit_matrix_path = get_file_path(
                dih_fit_matrix_file, dir_str=run_name, dir_path=data_path
            )

            dih_pred_matrix_path = get_file_path(
                dih_pred_matrix_file, dir_str=run_name, dir_path=data_path
            )

            rmsd_fit_matrix_path = get_file_path(
                rmsd_fit_matrix_file, dir_str=run_name, dir_path=data_path
            )

            rmsd_pred_matrix_path = get_file_path(
                rmsd_pred_matrix_file, dir_str=run_name, dir_path=data_path
            )

            nuc_df = mask_equal(df, nuc_class_col, nuc_class)

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
                included_df=dih_df,
                rmsd_matrix_path=rmsd_matrix_path,
                sup_resids=sup_resids,
                rmsd_resids=loop_resids,
                rmsd_atomids="CA",
                pair_aln=False,
                rmsd_dict=rmsd_dict,
                rmsd_json_path=rmsd_json_path,
            )

            build_dih_matrix(
                included_df=dih_df,
                max_norm_path=dih_matrix_path,
            )

            dih_matrix = load_matrix(dih_matrix_path)
            rmsd_matrix = load_matrix(rmsd_matrix_path)

            mask_dih_data(
                df=dih_df,
                matrix=dih_matrix,
                included_table_path=included_table_path,
                fit_matrix_path=dih_fit_matrix_path,
                removed_table_path=removed_table_path,
                pred_matrix_path=dih_pred_matrix_path,
                edia_dict=edia_dict,
                edia_min=0.4,
                edia_atomids="O",
            )
            mask_dih_data(
                df=dih_df,
                matrix=rmsd_matrix,
                included_table_path=included_table_path,
                fit_matrix_path=rmsd_fit_matrix_path,
                removed_table_path=removed_table_path,
                pred_matrix_path=rmsd_pred_matrix_path,
                edia_dict=edia_dict,
                edia_min=0.4,
                edia_atomids="O",
            )

            included_df = load_table(included_table_path)
            dih_fit_matrix = load_matrix(dih_fit_matrix_path)
            rmsd_fit_matrix = load_matrix(rmsd_fit_matrix_path)

            cluster_matrix(
                df=included_df,
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

            removed_df = load_table(removed_table_path)

            dih_pred_matrix = load_matrix(dih_pred_matrix_path)
            rmsd_pred_matrix = load_matrix(rmsd_pred_matrix_path)

            classify_matrix(
                cluster_df=cluster_df,
                removed_df=removed_df,
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
                    if rama in list(name_dict[loop_name][nuc_class].keys()):
                        name = name_dict[loop_name][nuc_class][rama]
                    else:
                        name = f"{loop_name}.{nuc_class}-Unknown-{i}"
                        i += 1

                rename_dict[cluster] = name

            result_df[cluster_col] = result_df[cluster_col].map(rename_dict)
            sum_df[cluster_col] = sum_df[cluster_col].map(rename_dict)
            classify_report_df[cluster_col] = cluster_report_df[cluster_col].map(
                rename_dict
            )

            result_df[loop_col] = loop_name
            sum_df[loop_col] = loop_name
            cluster_report_df[loop_col] = loop_name
            classify_report_df[loop_col] = loop_name

            result_df[nuc_class_col] = nuc_class
            sum_df[nuc_class_col] = nuc_class
            cluster_report_df[nuc_class_col] = nuc_class
            classify_report_df[nuc_class_col] = nuc_class

            loop_result_df = pd.concat([loop_result_df, result_df], sort=False)
            loop_sum_df = pd.concat([loop_result_df, sum_df], sort=False)
            loop_cluster_report_df = pd.concat(
                [loop_cluster_report_df, cluster_report_df], sort=False
            )
            loop_classify_report_df = pd.concat(
                [loop_classify_report_df, classify_report_df], sort=False
            )

        loop_result_df = loop_result_df.reset_index(drop=True)

        loop_sum_df = loop_sum_df.sort_values(
            by=[nuc_class_col, total_chain_col], ascending=[True, False]
        )

        loop_result_table_path = get_file_path(
            result_table_file, dir_str=loop_name, dir_path=data_path
        )
        loop_sum_table_path = get_file_path(
            sum_table_file, dir_str=loop_name, dir_path=data_path
        )
        loop_cluster_report_table_path = get_file_path(
            cluster_report_table_file, dir_str=loop_name, dir_path=data_path
        )
        loop_classify_report_table_path = get_file_path(
            classify_report_table_file, dir_str=loop_name, dir_path=data_path
        )

        save_table(loop_result_table_path, loop_result_df)
        save_table(loop_sum_table_path, loop_sum_df)
        save_table(loop_cluster_report_table_path, loop_cluster_report_df)
        save_table(loop_classify_report_table_path, loop_classify_report_df)

        cluster_dict = make_dict(
            lst_col(loop_result_df, pdb_id_col), lst_col(loop_result_df, cluster_col)
        )

        df[loop_name] = df[pdb_id_col].map(cluster_dict)

    save_table(entry_table_path, df)

    df = load_table(entry_table_path)

    gtp_atomids = ["O1G", "O2G", "O3G", "S1G", "S2G", "S3G"]

    dist_df = build_dist_table(
        mask_equal(df, sw1_name, sw1_gtp_name),
        x_resids=[32],
        y_resids=[bio_lig_col],
        x_atomids=["OH"],
        y_atomids=[gtp_atomids],
        atom_dist_col_lst=[atom_dist_col],
        hb_status_col_lst=[hb_status_col],
        outlier_col_lst=[outlier_col],
        check_hb=True,
    )

    sw1_gtp_dict = make_dict(
        lst_col(dist_df, pdb_id_col), lst_col(dist_df, hb_status_col)
    )

    df[hb_status_col] = df[pdb_id_col].map(sw1_gtp_dict).fillna("")
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

    print("Rascore clustering complete!")