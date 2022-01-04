# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from functions import *


def add_edia_status(
    df,
    edia_dict,
    edia_min=0.4,
    edia_resids=None,
    edia_atomids=None,
):
    if edia_atomids is None:
        edia_atomids = "O"

    edia_resid_lst = res_to_lst(edia_resids)
    edia_atomid_lst = type_lst(edia_atomids)

    for index in list(df.index.values):

        pdb_code = df.at[index, pdb_code_col]
        chainid = df.at[index, chainid_col]

        if edia_resid_lst is None:
            index_resid_lst = str_to_lst(df.at[index, bb_resid_col])
        else:
            index_resid_lst = list()
            for resid in edia_resid_lst:
                index_resid_lst.append(resid)

        edia_status = "Not Available"
        edia_below = False
        if pdb_code in list(edia_dict.keys()):
            if chainid in list(edia_dict[pdb_code].keys()):
                edia_status = f"Above {edia_min}"
                for index_resid in index_resid_lst:
                    if edia_below:
                        break
                    if index_resid in list(edia_dict[pdb_code][chainid].keys()):
                        for atomid in edia_atomid_lst:
                            if atomid in list(
                                edia_dict[pdb_code][chainid][index_resid].keys()
                            ):
                                edia_score = edia_dict[pdb_code][chainid][index_resid][
                                    atomid
                                ][edia_col]
                                if edia_score < edia_min:
                                    edia_below = True
                                    break

        if edia_below:
            edia_status = edia_status.replace("Above", "Below")

        df.at[index, edia_col] = edia_status

    return df


def mask_dih_table(df):

    included_df = mask_equal(df, complete_col, "True", reset_index=False)

    if edia_col in list(included_df.columns):
        included_df = mask_search(
            included_df, edia_col, "Below", " ", equal=False, reset_index=False
        )

    removed_df = df.loc[~df.index.isin(list(included_df.index.values)), :]

    return included_df, removed_df


def mask_dih_data(
    df,
    matrix,
    included_table_path,
    fit_matrix_path,
    removed_table_path=None,
    pred_matrix_path=None,
    edia_dict=None,
    edia_min=None,
    edia_resids=None,
    edia_atomids=None,
):

    if edia_dict is not None and pdb_code_col in list(df.columns):
        df = add_edia_status(
            df,
            edia_dict,
            edia_min=edia_min,
            edia_resids=edia_resids,
            edia_atomids=edia_atomids,
        )

    included_df, removed_df = mask_dih_table(df)

    included_df = order_rows(included_df, reset_index=False)
    removed_df = order_rows(removed_df, reset_index=False)

    included_index_lst = list(included_df.index.values)
    removed_index_lst = list(removed_df.index.values)

    included_df = included_df.reset_index(drop=True)
    save_table(included_table_path, included_df)

    fit_matrix = mask_matrix(matrix, included_index_lst, included_index_lst)
    save_matrix(fit_matrix_path, fit_matrix)

    if removed_table_path is not None:
        removed_df = removed_df.reset_index(drop=True)
        save_table(removed_table_path, removed_df)

    if pred_matrix_path is not None:
        pred_matrix = mask_matrix(matrix, removed_index_lst, included_index_lst)
        save_matrix(pred_matrix_path, pred_matrix)

    print("Masked dihedral data!")
