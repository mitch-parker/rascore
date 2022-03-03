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

import imp
import warnings

warnings.filterwarnings("ignore")

import pandas as pd
from tqdm import tqdm
import concurrent.futures

from ..functions.lst import lst_to_str, str_to_lst, res_to_lst
from ..functions.coord import build_add_resid_lst, resname_to_letter
from ..functions.col import (
    core_path_col,
    modelid_col,
    chainid_col,
    phi_col,
    psi_col,
    omega_col,
    rama_col,
    rotamer_col,
    sc_col_lst,
    complete_col,
    resname_col,
    dih_col_lst,
    bb_seq_col,
    bb_len_col,
    bb_resid_col,
    sc_len_col,
    sc_resid_col,
    sc_seq_col,
)
from ..functions.dih import get_rama_type, get_rot_type
from ..functions.table import (
    get_df_at_index,
    fix_val,
    get_val_col,
    get_col_val_lst,
    get_val_col_lst,
)
from ..functions.path import save_table


def get_index_dih(
    df, index, dih_dict, resid_dict, max_ca_dist=4, ext_mult=1, coord_path_col=None
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    index_df = get_df_at_index(df, index)

    coord_path = index_df.at[index, coord_path_col]

    if coord_path in list(dih_dict.keys()):

        modelid = fix_val(index_df.at[index, modelid_col], return_int=True)
        chainid = index_df.at[index, chainid_col]

        try:
            dih_resid_lst = list(dih_dict[coord_path][str(modelid)][chainid].keys())
        except:
            dih_resid_lst = list(dih_dict[coord_path][int(modelid)][chainid].keys())

        index_df.at[index, complete_col] = True

        bb_resid_lst = list()
        bb_seq_lst = list()

        sc_seq_lst = list()
        sc_resid_lst = list()
        sc_len_lst = list()

        for dih_col in dih_col_lst:

            resid_lst = resid_dict[dih_col]

            if resid_lst is not None:

                add_resid_lst = build_add_resid_lst(
                    coord_path,
                    modelid,
                    chainid,
                    resid_lst,
                    dih_resid_lst,
                    max_ca_dist=max_ca_dist,
                    ext_mult=ext_mult,
                )

                val = 999.00
                resname = "-"
                for i, add_resid in enumerate(add_resid_lst):

                    if str(add_resid) in dih_resid_lst:

                        try:
                            val = dih_dict[coord_path][str(modelid)][chainid][
                                str(add_resid)
                            ][dih_col]
                            resname = resname_to_letter(
                                dih_dict[coord_path][str(modelid)][chainid][
                                    str(add_resid)
                                ][resname_col]
                            )
                        except:
                            val = dih_dict[coord_path][int(modelid)][chainid][
                                str(add_resid)
                            ][dih_col]
                            resname = resname_to_letter(
                                dih_dict[coord_path][int(modelid)][chainid][
                                    str(add_resid)
                                ][resname_col]
                            )

                    if dih_col == phi_col:
                        bb_resid_lst.append(add_resid)
                        bb_seq_lst.append(resname)

                    if dih_col in sc_col_lst:
                        sc_resid_lst.append(add_resid)
                        sc_seq_lst.append(resname)

                    val_col = get_val_col(dih_col, i + 1)

                    index_df.at[index, val_col] = float(val)

                    if val == 999.00:
                        index_df.at[index, complete_col] = False

                    val = 999.00
                    resname = "-"

                if dih_col in sc_col_lst:
                    sc_len_lst.append(len(sc_resid_lst))

        index_df.at[index, bb_len_col] = len(bb_resid_lst)
        index_df.at[index, sc_len_col] = lst_to_str(sc_len_lst)
        index_df.at[index, bb_resid_col] = lst_to_str(bb_resid_lst)
        index_df.at[index, sc_resid_col] = lst_to_str(sc_resid_lst)
        index_df.at[index, bb_seq_col] = lst_to_str(bb_seq_lst, join_txt="")
        index_df.at[index, sc_seq_col] = lst_to_str(sc_seq_lst, join_txt="")

    else:
        index_df = pd.DataFrame()

    return index_df


def get_rama_str(bb_vals_lst):

    rama_str = ""

    for bb_vals in bb_vals_lst:

        if 999.00 in bb_vals:
            rama_type = "-"
        else:
            rama_type = get_rama_type(bb_vals[0], bb_vals[1], bb_vals[2])

        rama_str += rama_type

    return rama_str


def get_rot_str(rot_val_lst):

    rot_str = ""

    for rot_val in rot_val_lst:

        rot_str += get_rot_type(rot_val)

    return rot_str


def add_bb_rama(df):

    bb_resid_lst = get_col_val_lst(df, phi_col)

    if len(bb_resid_lst) > 0:

        index_lst = list(df.index.values)

        for index in index_lst:

            max_index = int(float(df.at[index, bb_len_col]))

            if max_index > 1:
                resid_lst = bb_resid_lst[:max_index]
            else:
                resid_lst = bb_resid_lst

            bb_val_lst = list()

            for resid in resid_lst[:max_index]:

                bb_val_lst.append(
                    (
                        df.at[index, get_val_col(phi_col, resid)],
                        df.at[index, get_val_col(psi_col, resid)],
                        df.at[index, get_val_col(omega_col, resid)],
                    )
                )

            df.at[index, rama_col] = get_rama_str(bb_val_lst)

    return df


def add_sc_rot(df):

    for index in list(df.index.values):

        if sc_len_col in list(df.columns):

            sc_len_lst = str_to_lst(df.at[index, sc_len_col])

            sc_val_lst = list()
            i = 0

            for sc_col in sc_col_lst:

                sc_resid_lst = get_col_val_lst(df, sc_col)

                if len(sc_resid_lst) > 0:

                    max_index = int(float(sc_len_lst[i]))

                    if max_index > 1:
                        resid_lst = sc_resid_lst[:max_index]
                    else:
                        resid_lst = sc_resid_lst

                    for resid in resid_lst:

                        sc_val_lst.append(df.at[index, get_val_col(sc_col, resid)])

                    i += 1

            if len(sc_val_lst) > 0:

                df.at[index, rotamer_col] = get_rot_str(sc_val_lst)

    return df


def build_dih_table(
    df,
    dih_dict,
    dih_table_path=None,
    bb_resids=None,
    chi1_resids=None,
    chi2_resids=None,
    altchi1_resids=None,
    altchi2_resids=None,
    chi3_resids=None,
    chi4_resids=None,
    chi5_resids=None,
    max_ca_dist=4,
    ext_mult=1,
    coord_path_col=None,
    num_cpu=1,
):
    df = df.reset_index(drop=True)

    df_col_lst = list(df.columns)

    check_col_lst = [
        bb_resid_col,
        sc_resid_col,
        bb_seq_col,
        sc_seq_col,
        bb_len_col,
        sc_len_col,
        rama_col,
        rotamer_col,
    ]

    for check_col in check_col_lst:
        if check_col in df_col_lst:
            del df[check_col]

    for dih_col in dih_col_lst:
        for val_col in get_val_col_lst(df, dih_col):
            del df[val_col]

    angle_resids_lst = [
        bb_resids,
        bb_resids,
        bb_resids,
        chi1_resids,
        chi2_resids,
        altchi1_resids,
        altchi2_resids,
        chi3_resids,
        chi4_resids,
        chi5_resids,
    ]

    resid_dict = dict()

    for index, angle_resids in enumerate(angle_resids_lst):

        resid_dict[dih_col_lst[index]] = res_to_lst(angle_resids)

    dih_df = pd.DataFrame()

    if num_cpu == 1:
        for index in tqdm(
            list(df.index.values),
            desc="Building dihedral table",
            position=0,
            leave=True,
        ):
            index_df = get_index_dih(
                df,
                index,
                dih_dict,
                resid_dict,
                max_ca_dist=max_ca_dist,
                ext_mult=ext_mult,
                coord_path_col=coord_path_col,
            )
            for col in list(index_df.columns):
                if col not in list(dih_df.columns):
                    dih_df[col] = 999.00

            dih_df = pd.concat([dih_df, index_df], sort=False)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    get_index_dih,
                    df,
                    index,
                    dih_dict,
                    resid_dict,
                    max_ca_dist=max_ca_dist,
                    ext_mult=ext_mult,
                    coord_path_col=coord_path_col,
                )
                for index in list(df.index.values)
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Building dihedral table",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):
                index_df = job.result()

                for col in list(index_df.columns):
                    if col not in list(dih_df.columns):
                        dih_df[col] = 999.00

                dih_df = pd.concat([dih_df, index_df], sort=False)

    dih_df = dih_df.reset_index(drop=True)

    for dih_col in dih_col_lst:
        for val_col in get_val_col_lst(dih_df, dih_col):
            dih_df[val_col] = dih_df[val_col].fillna(999.00)

    dih_df = add_bb_rama(dih_df)
    dih_df = add_sc_rot(dih_df)

    print("Built dihedral table!")

    if dih_table_path is not None:
        save_table(dih_table_path, dih_df)
    else:
        return dih_df