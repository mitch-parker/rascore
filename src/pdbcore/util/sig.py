# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import pandas as pd
import numpy as np

import statsmodels.api as sm
from statsmodels.sandbox.stats.multicomp import multipletests

from util.lst import type_lst
from util.data import lst_col, mask_equal, mask_unequal
from util.col import (
    total_col,
    index_col,
    p_col,
    correct_p_col,
    a_col,
    b_col,
    c_col,
    d_col,
    risk_ratio_col,
    low_ci_col,
    up_ci_col,
    sig_col,
)


def label_sig(val, label="*", cutoff=0.05):

    if val > cutoff:
        return "ns"
    else:
        return label


def calc_rr(df, exp_cols, out_cols, correct_method="fdr_bh"):

    df = df.dropna()
    df = df.reset_index(drop=True)

    i_col_lst = type_lst(out_cols)
    j_col_lst = type_lst(exp_cols)

    col_lst = i_col_lst + j_col_lst

    if total_col not in list(df.columns):

        df.index.name = index_col
        df = df.reset_index()

        df = df.groupby(col_lst)[index_col].nunique().reset_index(name=total_col)

    total = df[total_col].sum()

    df_index_lst = list(df.index.values)

    for index in df_index_lst:

        i_df = df.copy(deep=True)
        j_df = df.copy(deep=True)

        for col in i_col_lst:

            val = df.at[index, col]

            i_df = mask_equal(i_df, col, val)

        for col in j_col_lst:

            val = df.at[index, col]

            j_df = mask_equal(j_df, col, val)

        a = df.at[index, total_col]
        b = i_df[total_col].sum() - a
        c = j_df[total_col].sum() - a
        d = total - (a + b + c)

        matrix = np.zeros((2, 2))

        matrix[0, 0] = a
        matrix[0, 1] = b
        matrix[1, 0] = c
        matrix[1, 1] = d

        table = sm.stats.Table2x2(matrix)

        risk_ratio = table.riskratio
        p_val = table.riskratio_pvalue()
        low_ci = table.riskratio_confint()[0]
        up_ci = table.riskratio_confint()[1]

        df.at[index, a_col] = a
        df.at[index, b_col] = b
        df.at[index, c_col] = c
        df.at[index, d_col] = d
        df.at[index, risk_ratio_col] = risk_ratio
        df.at[index, p_col] = p_val
        df.at[index, low_ci_col] = low_ci
        df.at[index, up_ci_col] = up_ci

    if correct_method is not None:

        p_lst = lst_col(df, p_col)

        correct_p_lst = multipletests(p_lst, method=correct_method)[1]

        for index in df_index_lst:

            correct_p = correct_p_lst[index]

            df.at[index, correct_p_col] = correct_p

    df = df.sort_values(by=col_lst)

    noise_df = pd.DataFrame()
    for col in col_lst:
        temp_df = mask_equal(df, col, "Noise")
        noise_df = pd.concat([noise_df, temp_df], sort=False)
        df = mask_unequal(df, col, "Noise")

    df = pd.concat([df, noise_df], sort=False)

    df = df.reset_index(drop=True)

    df[sig_col] = df[correct_p_col].map(label_sig)

    return df