# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from .scripts import *
from .constants import *


def prep_table(file_name, data_path=None):

    file_path = get_file_path(file_name, dir_path=data_path)

    df = load_table(file_path)

    if df is not None:
        df_col_lst = list(df.columns)
        for index in list(df.index.values):
            for col in [x for x in path_col_lst if x in df_col_lst]:
                dir_str = col.split("_path")[0]
                df.at[index, col] = get_file_path(
                    get_file_name(df.at[index, col]), dir_path=f"{data_path}/{dir_str}"
                )

        save_table(file_path, df)


def prep_json(file_name, data_path=None):

    file_path = get_file_path(file_name, dir_path=data_path)

    old_dict = load_json(file_path)

    if old_dict is not None:
        new_dict = dict()
        for old_key in list(old_dict.keys()):
            for col in path_col_lst:
                dir_str = col.split("_path")[0]
                if dir_str in old_key:
                    break
            new_key = get_file_path(
                get_file_name(old_key), dir_path=f"{data_path}/{dir_str}"
            )
            new_dict[new_key] = old_dict[old_key]
            for sub_key in list(new_dict[new_key].keys()):
                for col in path_col_lst:
                    if col in list(new_dict[new_key][sub_key].keys()):
                        new_dict[new_key][sub_key][col] = get_file_path(
                            get_file_name(new_dict[new_key][sub_key][col]),
                            dir_path=f"{data_path}/{dir_str}",
                        )

        save_json(file_path, new_dict)


def prep_rascore(data_path=None):

    table_file_lst = [entry_table_file, interf_table_file, pocket_table_file]

    json_file_lst = [interf_json_file, pocket_json_file, dih_json_file]

    for table_file in tqdm(
        table_file_lst,
        desc="Preparing rascore database tables",
        position=0,
        leave=True,
    ):
        prep_table(table_file, data_path=data_path)

    for json_file in tqdm(
        json_file_lst,
        desc="Preparing rascore database jsons",
        position=0,
        leave=True,
    ):
        prep_json(json_file, data_path=data_path)
