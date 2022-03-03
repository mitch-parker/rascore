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

import os
from tqdm import tqdm

from ..functions.path import (
    get_file_path,
    load_table,
    get_file_name,
    save_table,
    save_json,
    load_json,
    rascore_str,
    build_str,
)
from ..functions.file import (
    entry_table_file,
    interf_table_file,
    pocket_table_file,
    interf_json_file,
    pocket_json_file,
    dih_json_file,
)
from ..functions.col import path_col_lst


def prep_table(file_name, build_path=None):

    file_path = get_file_path(file_name, dir_path=build_path)

    df = load_table(file_path)

    if df is not None:
        df_col_lst = list(df.columns)
        for index in list(df.index.values):
            for col in [x for x in path_col_lst if x in df_col_lst]:
                dir_str = col.split("_path")[0]
                df.at[index, col] = get_file_path(
                    get_file_name(df.at[index, col]), dir_path=f"{build_path}/{dir_str}"
                )

        save_table(file_path, df)


def prep_json(file_name, build_path=None):

    file_path = get_file_path(file_name, dir_path=build_path)

    old_dict = load_json(file_path)

    if old_dict is not None:
        new_dict = dict()
        for old_key in list(old_dict.keys()):
            for col in path_col_lst:
                dir_str = col.split("_path")[0]
                if f"/{dir_str}/" in old_key:
                    break
            new_key = get_file_path(
                get_file_name(old_key), dir_path=f"{build_path}/{dir_str}"
            )
            new_dict[new_key] = old_dict[old_key]
            for sub_key in list(new_dict[new_key].keys()):
                for col in path_col_lst:
                    if col in list(new_dict[new_key][sub_key].keys()):
                        sub_dir_str = col.split("_path")[0]
                        new_dict[new_key][sub_key][col] = get_file_path(
                            get_file_name(new_dict[new_key][sub_key][col]),
                            dir_path=f"{build_path}/{sub_dir_str}",
                        )

        save_json(file_path, new_dict)


def prep_rascore(build_path=None):

    if build_path is None:
        build_path = f"{os.getcwd()}/{rascore_str}_{build_str}"

    table_file_lst = [entry_table_file, interf_table_file, pocket_table_file]

    json_file_lst = [interf_json_file, pocket_json_file, dih_json_file]

    for table_file in tqdm(
        table_file_lst,
        desc="Preparing Rascore database tables",
        position=0,
        leave=True,
    ):
        prep_table(table_file, build_path=build_path)

    for json_file in tqdm(
        json_file_lst,
        desc="Preparing Rascore database jsons",
        position=0,
        leave=True,
    ):
        prep_json(json_file, build_path=build_path)
