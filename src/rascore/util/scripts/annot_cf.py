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
from tqdm import tqdm
import concurrent.futures

from ..functions.lst import type_lst
from ..functions.coord import load_cif_dict, search_cif_dict
from ..functions.path import save_table
from ..functions.table import merge_tables
from ..functions.col import (
    len_a_col,
    len_b_col,
    len_c_col,
    ang_a_col,
    ang_b_col,
    ang_g_col,
    rcsb_path_col,
    space_col,
    cf_col,
)

len_col_lst = [len_a_col, len_b_col, len_c_col]
ang_col_lst = [ang_a_col, ang_b_col, ang_g_col]
len_letter_lst = ["a", "b", "c"]
ang_letter_lst = ["alpha", "beta", "gamma"]


def get_path_len_ang_df(coord_path):

    df = pd.DataFrame()

    cif_dict = load_cif_dict(coord_path)

    df.at[0, rcsb_path_col] = coord_path

    if search_cif_dict(cif_dict, "_exptl.method") == "X-RAY DIFFRACTION":

        df.at[0, space_col] = search_cif_dict(
            cif_dict, "_symmetry.space_group_name_H-M"
        )

        for col, letter in zip(len_col_lst, len_letter_lst):

            df.at[0, col] = search_cif_dict(cif_dict, f"_cell.length_{letter}")

        for col, letter in zip(ang_col_lst, ang_letter_lst):

            df.at[0, col] = search_cif_dict(cif_dict, f"_cell.angle_{letter}")

    else:
        df.at[0, space_col] = "None"
        for col, letter in zip(len_col_lst, len_letter_lst):
            df.at[0, col] = 999.00
        for col, letter in zip(ang_col_lst, ang_letter_lst):
            df.at[0, col] = 999.00

    return df


def add_cf(df, min_simi=0.05):

    val_dict = dict()
    prev_dict = dict()

    col_lst = len_col_lst + ang_col_lst

    for index in tqdm(
        list(df.index.values), desc="Annotating crystal forms", position=0, leave=True
    ):

        curr_path = df.at[index, rcsb_path_col]
        curr_space = df.at[index, space_col]

        if curr_space == "None":
            df.at[index, cf_col] = "None"
        else:
            val_dict[curr_path] = {
                len_a_col: df.at[index, len_a_col],
                len_b_col: df.at[index, len_b_col],
                len_c_col: df.at[index, len_c_col],
                ang_a_col: df.at[index, ang_a_col],
                ang_b_col: df.at[index, ang_b_col],
                ang_g_col: df.at[index, ang_g_col],
            }

            if curr_space not in list(prev_dict.keys()):
                id = 1
                prev_dict[curr_space] = {id: type_lst(curr_path)}

            else:

                new_form = True

                for id in list(prev_dict[curr_space].keys()):

                    prev_path_lst = prev_dict[curr_space][id]

                    diff_form = False

                    for col in col_lst:

                        diff_lst = list()

                        for prev_path in prev_path_lst:

                            curr_val = float(val_dict[curr_path][col])
                            prev_val = float(val_dict[prev_path][col])

                            diff_lst.append(abs(prev_val - curr_val) / prev_val)

                        mean_diff = np.mean(diff_lst)

                        if mean_diff > min_simi:
                            diff_form = True
                            break

                    if not diff_form:

                        prev_dict[curr_space][id].append(curr_path)

                        new_form = False

                        break

                if new_form:

                    id = max(list(prev_dict[curr_space].keys())) + 1

                    prev_dict[curr_space][id] = type_lst(curr_path)

            df.at[index, cf_col] = f"{curr_space} (CF{id})"

    return df


def annot_cf(coord_paths, cf_table_path=None, min_simi=0.05, data=None, num_cpu=1):

    coord_path_lst = type_lst(coord_paths)

    df = pd.DataFrame()

    if num_cpu == 1:
        for coord_path in tqdm(
            coord_path_lst, desc="Getting crystal information", position=0, leave=True
        ):

            df = pd.concat([df, get_path_len_ang_df(coord_path)], sort=False)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(get_path_len_ang_df, coord_path)
                for coord_path in coord_path_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Getting crystal information",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                df = pd.concat([df, job.result()], sort=False)

    df = df.reset_index(drop=True)

    df = add_cf(df, min_simi=min_simi)

    if data is not None:
        df_col_lst = list(data.columns)

        for col in [cf_col, space_col, len_a_col, len_b_col, len_c_col, ang_a_col, ang_b_col, ang_g_col]:
            if col in df_col_lst:
                del data[col]

        df = merge_tables(df, data)

    print("Annotated crystal forms!")

    if cf_table_path is not None:
        save_table(cf_table_path, df)
    else:
        return df