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

import re
import pandas as pd
import numpy as np

from .col import (
    reformat_col_lst,
    resid_col,
    order_col_lst,
    interf_path_col,
    pocket_path_col,
    core_path_col,
    pdb_id_col,
    pdb_code_col,
    modelid_col,
    chainid_col,
    interf_col,
    cf_col,
    index_col,
    total_col,
    prot_col,
    count_col_dict,
    data_col_lst,
)
from .lst import (
    format_val,
    str_to_lst,
    type_lst,
    sort_lst,
    move_end_lst,
    lst_nums,
    lst_unique,
    lst_to_str,
    lst_inter,
)


def fix_val(val, return_str=False, return_int=False, return_float=False):

    if type(val) == str and is_int(val):
        if val[-2:] == ".0":
            val = str(val).replace(".0", "")

    return format_val(
        val, return_str=return_str, return_int=return_int, return_float=return_float
    )


def fix_query(query):

    return [fix_val(x) for x in type_lst(query)]


def fix_col(df, col):

    for index in list(df.index.values):

        df.at[index, col] = fix_val(df.at[index, col])

    return df


def mask_equal(df, col, query, reset_index=True):

    temp_df = df.copy(deep=True)

    query = fix_query(query)

    temp_df = fix_col(temp_df, col)

    temp_df = temp_df[temp_df[col].isin(query)]

    if reset_index:
        temp_df = temp_df.reset_index(drop=True)

    return temp_df


def mask_unequal(df, col, query, reset_index=True):

    temp_df = df.copy(deep=True)

    query = fix_query(query)

    temp_df = fix_col(temp_df, col)

    temp_df = temp_df[~temp_df[col].isin(query)]

    if reset_index:
        temp_df = temp_df.reset_index(drop=True)

    return temp_df


def mask_search(df, col, items, sep_txt, equal=True, reset_index=True):

    temp_df = df.copy(deep=True)

    item_lst = fix_query(items)

    temp_df = fix_col(temp_df, col)

    bool_lst = list()

    for row in df[col]:
        row_lst = str_to_lst(row, sep_txt=sep_txt, return_str=True)
        in_lst = not equal
        for item in item_lst:
            if str(item) in row_lst:
                in_lst = equal
        bool_lst.append(in_lst)

    temp_df = temp_df[bool_lst]

    if reset_index:
        temp_df = temp_df.reset_index(drop=True)

    return temp_df


def mask_greater(df, col, cutoff, reset_index=True):

    temp_df = df.copy(deep=True)

    temp_df[col] = temp_df[col].map(float)

    mask = temp_df[col] >= cutoff
    temp_df = temp_df.loc[mask, :]

    if reset_index:
        temp_df = temp_df.reset_index(drop=True)

    return temp_df


def mask_less(df, col, cutoff, reset_index=True):

    temp_df = df.copy(deep=True)

    temp_df[col] = temp_df[col].map(float)

    mask = temp_df[col] <= cutoff
    temp_df = temp_df.loc[mask, :]

    if reset_index:
        temp_df = temp_df.reset_index(drop=True)

    return temp_df


def mask_between(df, col, bottom, top, reset_index=True):

    return mask_greater(
        mask_less(df, col, top, reset_index=reset_index),
        col,
        bottom,
        reset_index=reset_index,
    )


def mask_matrix(matrix, row_lst, col_lst):

    return matrix[np.ix_(row_lst, col_lst)]


def lst_col(
    df, col, unique=False, return_str=False, return_int=False, return_float=False
):

    if unique:
        val_lst = list(df[col].unique())
        val_lst = sort_lst(val_lst)
        val_lst = move_end_lst(val_lst, ["Noise", "None"])
    else:
        val_lst = list(df[col].to_list())

    val_lst = format_val(
        val_lst,
        return_str=return_str,
        return_int=return_int,
        return_float=return_float,
    )

    return val_lst


def get_val_col(col, resid):

    return f"{col}_{resid}"


def get_col_val(col):

    return int(col.rsplit("_", 1)[1])


def get_col_col_lst(df, col):

    return [x for x in list(df.columns) if col in x]


def get_col_val_lst(df, col):

    return [get_col_val(x) for x in get_col_col_lst(df, col)]


def get_val_col_lst(df, val):

    return [x for x in list(df.columns) if val in x]


def reformat_val_table(df, cols):

    col_lst = type_lst(cols)

    id_lst = list()

    for col in list(df.columns):
        add = True
        for reformat_col in reformat_col_lst:
            if reformat_col in col:
                add = False
                break
        if add:
            id_lst.append(col)

    temp_df_lst = list()

    for col in col_lst:

        val_col_lst = get_val_col_lst(df, col)
        col_val_lst = get_col_val_lst(df, col)

        col_val_dict = make_dict(val_col_lst, col_val_lst)

        temp_df = pd.melt(
            df,
            id_vars=id_lst,
            value_vars=val_col_lst,
            var_name=resid_col,
            value_name=col,
        ).drop_duplicates()

        temp_df[resid_col] = temp_df[resid_col].map(col_val_dict)

        temp_df_lst.append(temp_df)

    if len(temp_df_lst) == 1:
        df = temp_df_lst[0]
    else:
        df = merge_tables(temp_df_lst[0], temp_df_lst[1])

        if len(temp_df_lst) > 2:
            i_lst = lst_nums(2, len(temp_df_lst) - 1)

            for i in i_lst:
                df = merge_tables(df, temp_df_lst[i])

    return df


def get_col_order(df):

    df_col_lst = list(df.columns)

    col_lst = list()
    for order_col in order_col_lst:
        if order_col in df_col_lst:
            col_lst.append(order_col)
        if order_col in data_col_lst:
            for col in get_val_col_lst(df, order_col):
                if col not in col_lst:
                    col_lst.append(col)

    for col in df_col_lst:
        if col not in col_lst:
            col_lst.append(col)

    return col_lst


def order_cols(df, col_lst):

    df = df.reindex(columns=col_lst)

    return df


def order_rows(df, col_lst=None, reset_index=False):

    df_col_lst = list(df.columns)

    if col_lst is None:
        col_lst = list()

        if interf_path_col in df_col_lst:
            col_lst.append(interf_path_col)
        elif pocket_path_col in df_col_lst:
            col_lst.append(pocket_path_col)
        elif core_path_col in df_col_lst:
            col_lst.append(core_path_col)

        if pdb_id_col in df_col_lst:
            col_lst.append(pdb_id_col)
        if pdb_code_col in df_col_lst:
            col_lst.append(pdb_code_col)
        if modelid_col in df_col_lst:
            col_lst.append(modelid_col)
        if chainid_col in df_col_lst:
            col_lst.append(chainid_col)
        if interf_col in df_col_lst:
            col_lst.append(interf_col)

    if len(col_lst) > 0:
        df = df.sort_values(by=col_lst)

    if reset_index:
        df = df.reset_index(drop=True)

    return df


def make_dict(lst_1, lst_2):

    lst_dict = dict(zip(lst_1, lst_2))

    return lst_dict


def merge_dicts(dict_lst):

    fin_dict = {}

    for curr_dict in dict_lst:

        fin_dict = {**fin_dict, **curr_dict}

    return fin_dict


def rename_dict_key(the_dict, old_key, new_key):

    the_dict[new_key] = the_dict.pop(old_key)

    return the_dict


def get_str_num(str_num):

    return int("".join(filter(str.isdigit, str_num)))


def replace_str(val, term, replace=""):

    term_lst = type_lst(term)

    for term in term_lst:
        if term in val:
            val = val.replace(term, replace)

    return val


def build_count_dict(val_lst):

    count_dict = {}
    for val in val_lst:
        if val in count_dict:
            count_dict[val] += 1
        else:
            count_dict[val] = 1

    return count_dict


def build_col_count_dict(df, col, col_lst=None, return_str=False, return_int=False):

    if col_lst is not None:
        col_lst = type_lst(col_lst)
        col_lst.append(col)
        df = df.loc[:, col_lst]
        df = df.drop_duplicates()

    count_lst = lst_col(df, col, return_str=return_str, return_int=return_int)
    count_dict = build_count_dict(count_lst)

    return count_dict


def lst_by_freq(val_lst):

    count_dict = build_count_dict(val_lst)

    sorted_vals = sorted(count_dict, key=count_dict.get, reverse=True)

    return sorted_vals


def get_col_most_common(df, col, n=None):

    col_lst = lst_unique(lst_by_freq(lst_col(df, col)))

    if n is not None:
        col_lst[:n]
        if n == 1:
            col_lst = col_lst[0]

    return col_lst


def extract_str(val):

    if val is not None:
        val = lst_to_str([x for x in str(val) if x.isalpha()], join_txt="", empty="")

    return val


def extract_int(val):

    if val is not None:
        val = int(re.search(r"\d+", str(val)).group())

    return val


def is_str(val):

    if type(val) == str:
        status = len(extract_str(val)) == len(val)
    else:
        status = False

    return status


def is_int(val):

    return not is_str(val)


def explode_table(df, cols, sep=","):

    temp_df = df.copy(deep=True)

    temp_df[cols] = temp_df[cols].str.split(sep)

    temp_df = temp_df.explode(cols)
    temp_df = temp_df.reset_index(drop=True)

    return temp_df


def merge_tables(left_df, right_df, how="left"):

    col_lst = lst_inter(list(left_df.columns), list(right_df.columns))

    for col in col_lst:
        left_df[col] = left_df[col].map(str)
        right_df[col] = right_df[col].map(str)

    return left_df.merge(
        right_df,
        on=col_lst,
        how=how,
    )


def build_label_dict(
    df,
    col,
    return_str=False,
    return_int=False,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
):

    df_col_lst = list(df.columns)

    if interf_path_col in df_col_lst:
        interf_dict = build_col_count_dict(
            df,
            col,
            col_lst=[interf_path_col],
            return_str=return_str,
            return_int=return_int,
        )
    if pocket_path_col in df_col_lst:
        pocket_dict = build_col_count_dict(
            df,
            col,
            col_lst=[pocket_path_col],
            return_str=return_str,
            return_int=return_int,
        )

    if count_chain:
        if pdb_id_col in df_col_lst:
            chain_col_lst = [pdb_id_col]
        else:
            chain_col_lst = [core_path_col, modelid_col, chainid_col]

        chain_dict = build_col_count_dict(
            df, col, col_lst=chain_col_lst, return_str=return_str, return_int=return_int
        )

    if count_pdb:
        if pdb_code_col in df_col_lst:
            entry_dict = build_col_count_dict(
                df,
                col,
                col_lst=[pdb_code_col],
                return_str=return_str,
                return_int=return_int,
            )

    if count_cf:
        if cf_col in df_col_lst:
            cf_dict = build_col_count_dict(
                df, col, col_lst=[cf_col], return_str=return_str, return_int=return_int
            )

    val_lst = lst_col(
        df, col, unique=True, return_str=return_str, return_int=return_int
    )

    label_dict = dict()

    for val in val_lst:
        name = val
        name += " ("

        if interf_path_col in df_col_lst:
            name += f"I={interf_dict[str(val)]}; "
        if pocket_path_col in df_col_lst:
            name += f"S={pocket_dict[str(val)]}; "

        if count_chain:
            name += f"N={chain_dict[str(val)]}"
        if count_pdb:
            if pdb_code_col in df_col_lst:
                if count_chain:
                    name += "; "
                name += f"PDB={entry_dict[str(val)]}"
        if count_cf:
            if cf_col in df_col_lst:
                if count_chain or (count_pdb and pdb_code_col in df_col_lst):
                    name += "; "
                name += f"CF={cf_dict[str(val)]}"

        name += ")"
        label_dict[val] = name

    return label_dict


def build_label_color_dict(
    df,
    col,
    color_dict,
    label_order=None,
    return_str=False,
    return_int=False,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
):

    label_dict = build_label_dict(
        df,
        col,
        return_str=return_str,
        return_int=return_int,
        count_chain=count_chain,
        count_pdb=count_pdb,
        count_cf=count_cf,
    )

    if label_order is None:
        label_order = list(label_dict.keys())

    label_color_dict = dict()
    for label in label_order:
        label_color_dict[label_dict[label]] = color_dict[label]

    return label_color_dict


def get_ncols(label_lst):

    total = len(label_lst)

    if 1 < total <= 5:
        cols = total
    else:
        cols = 5

    return cols


def title_str(str_val):

    if "_" in str_val:
        str_val = str_val.replace("_", " ")

    return str_val.title()


def get_val_index_lst(df, col, val):

    temp_df = mask_equal(df, col, val, reset_index=False)

    return list(temp_df.index.values)


def build_count_table(df, cols):

    col_lst = type_lst(cols)

    temp_df = df.copy(deep=True)

    temp_df = temp_df.rename(columns=count_col_dict)

    df_col_lst = list(temp_df.columns)

    count_df_lst = list()

    for col in list(count_col_dict.values()):
        if col in df_col_lst:
            count_df_lst.append(
                pd.pivot_table(temp_df, index=col_lst, values=col, aggfunc="nunique")
                .fillna("-")
                .reset_index()
            )
    total_counts = len(count_df_lst)
    if total_counts == 0:
        temp_df = temp_df.reset_index()
        temp_df = temp_df.rename(columns={index_col: total_col})
        count_df = (
            pd.pivot_table(temp_df, index=col_lst, values=total_col, aggfunc="nunique")
            .fillna("-")
            .reset_index()
        )
    else:
        count_df = count_df_lst[0]
        if total_counts > 0:
            for i in lst_nums(1, total_counts - 1):
                count_df = merge_tables(count_df, count_df_lst[i])

        if prot_col in df_col_lst:
            prot_df = pd.crosstab(
                columns=temp_df[prot_col], index=[temp_df[col] for col in col_lst]
            )
            for col in list(prot_df.columns):
                prot_df[col] = prot_df[col].replace({0: "-"})

            prot_df = prot_df.reset_index()

            count_df = merge_tables(count_df, prot_df)

    for col in col_lst:
        count_df[col] = count_df[col].map(str)

    return count_df


def get_df_at_index(df, index):

    return df.loc[type_lst(index), :]


def convert_col_percent(df, col, dec=1, return_frac=False):

    df[col] = (df[col] / df[col].sum())

    if not return_frac:
        df[col] = df[col] * 100


    df[col] = df[col].apply(lambda x: round(x, dec))

    return df

def str_to_dict(dict_str, return_str=False, return_int=False, return_float=False):

    str_dict = dict()

    for key_val in str_to_lst(dict_str ,sep_txt="|"):

        str_dict[key_val.split(":")[0]] = str_to_lst(key_val.split(":")[1], return_str=return_str, return_int=return_int, return_float=return_float)

    return str_dict
