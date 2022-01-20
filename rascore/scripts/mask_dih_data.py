# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from ..functions import *


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

    fit_df = mask_equal(df, complete_col, "True", reset_index=False)

    if edia_col in list(fit_df.columns):
        fit_df = mask_search(
            fit_df, edia_col, "Below", " ", equal=False, reset_index=False
        )

    pred_df = df.loc[~df.index.isin(list(fit_df.index.values)), :]

    return fit_df, pred_df


def mask_dih_data(
    df,
    matrix,
    fit_table_path,
    fit_matrix_path,
    pred_table_path=None,
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

    fit_df, pred_df = mask_dih_table(df)

    fit_df = order_rows(fit_df, reset_index=False)
    pred_df = order_rows(pred_df, reset_index=False)

    fit_index_lst = list(fit_df.index.values)
    pred_index_lst = list(pred_df.index.values)

    fit_df = fit_df.reset_index(drop=True)
    save_table(fit_table_path, fit_df)

    fit_matrix = mask_matrix(matrix, fit_index_lst, fit_index_lst)
    save_matrix(fit_matrix_path, fit_matrix)

    if pred_table_path is not None:
        pred_df = pred_df.reset_index(drop=True)
        save_table(pred_table_path, pred_df)

    if pred_matrix_path is not None:
        pred_matrix = mask_matrix(matrix, pred_index_lst, fit_index_lst)
        save_matrix(pred_matrix_path, pred_matrix)

    print("Masked dihedral data!")
