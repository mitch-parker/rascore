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
import numpy as np

import statsmodels.api as sm
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats

from .lst import type_lst
from .table import lst_col, mask_equal, mask_unequal
from .col import (
    p_col,
    correct_p_col,
    total_col,
    index_col,
    a_col,
    b_col,
    c_col,
    d_col,
    risk_ratio_col,
    up_ci_col,
    low_ci_col,
    sig_col,
    corr_col,
)


def label_sig(val, label="*", cutoff=0.05):

    if val > cutoff:
        return "ns"
    else:
        return label


def correct_p_vals(df, correct_method="fdr_bh"):

    df = df.reset_index(drop=True)

    p_lst = lst_col(df, p_col)

    correct_p_lst = multipletests(p_lst, method=correct_method)[1]

    for index in list(df.index.values):

        correct_p = correct_p_lst[index]

        df.at[index, correct_p_col] = correct_p

    return df


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
        df = correct_p_vals(df, correct_method=correct_method)

    df = df.sort_values(by=col_lst)

    noise_df = pd.DataFrame()
    for col in col_lst:
        temcorr_df = mask_equal(df, col, "Noise")
        noise_df = pd.concat([noise_df, temcorr_df], sort=False)
        df = mask_unequal(df, col, "Noise")

    df = pd.concat([df, noise_df], sort=False)

    df = df.reset_index(drop=True)

    df[sig_col] = df[correct_p_col].map(label_sig)

    return df


def calc_corr_stat(df, x_col, y_col, return_df=False, use_kt=False):

    if len(df) < 2:
        corr = np.nan
        p_val = np.nan
    else:
        df[x_col] = df[x_col].map(float)
        df[y_col] = df[y_col].map(float)

        x_lst = lst_col(df, x_col)
        y_lst = lst_col(df, y_col)

        if use_kt:
            corr, p_val = stats.kendalltau(x_lst, y_lst)
        else:
            corr, p_val = stats.pearsonr(x_lst, y_lst)

    if return_df:
        return pd.DataFrame({corr_col: [corr], p_col: [p_val]})
    else:
        return (corr, p_val)


def calc_corr(df, x_col, y_col, hue_cols=None, correct_method="fdr_bh", use_kt=False):

    val_df = df.copy(deep=True)

    if hue_cols is None:
        corr_df = calc_corr_stat(val_df, x_col, y_col, return_df=True, use_kt=use_kt)
    else:
        hue_col_lst = type_lst(hue_cols)

        corr_df = val_df.copy(deep=True)

        corr_df.index.name = index_col
        corr_df = corr_df.reset_index()

        corr_df = (
            corr_df.groupby(hue_col_lst)[index_col]
            .nunique()
            .reset_index(name=total_col)
        )

        df_index_lst = list(corr_df.index.values)

        for index in df_index_lst:

            hue_df = val_df.copy(deep=True)

            for col in hue_col_lst:

                val = corr_df.at[index, col]

                hue_df = mask_equal(hue_df, col, val)

            corr = calc_corr_stat(hue_df, x_col, y_col, use_kt=use_kt)

            corr_df.at[index, corr_col] = corr[0]
            corr_df.at[index, p_col] = corr[1]

        if correct_method is not None:
            corr_df = correct_p_vals(corr_df, correct_method=correct_method)

    corr_df = corr_df.sort_values(by=hue_col_lst)

    corr_df = corr_df.reset_index(drop=True)

    corr_df[sig_col] = corr_df[correct_p_col].map(label_sig)

    return corr_df
