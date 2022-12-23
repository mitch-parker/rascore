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
from tqdm import tqdm

from ..functions.coord import get_pdb_id
from ..functions.pdbaa import (
    get_pdbaa_swiss_id,
    get_pdbaa_chainid,
    get_pdbaa_method,
    get_pdbaa_pdb_code,
    get_pdbaa_prot,
    get_pdbaa_r_factor,
    get_pdbaa_resolution,
)
from ..functions.table import lst_col
from ..functions.seq import load_record_dict, get_record_desc, get_record_seq
from ..functions.col import (
    pdb_id_col,
    pdb_code_col,
    chainid_col,
    swiss_id_col,
    prot_col,
    method_col,
    resolution_col,
    r_factor_col,
    seq_col,
    len_col,
    swiss_id_col,
)
from ..functions.path import save_table


def search_pdbaa(
    pdbaa_fasta_path,
    search_lst,
    fix_dict=None,
    entry_table_path=None,
    min_length=None,
    min_resolution=None,
):

    pdbaa_dict = load_record_dict(pdbaa_fasta_path)

    for search in search_lst:
        if "_" in search:
            search_col = swiss_id_col
            break
        else:
            search_col = pdb_id_col
            break

    col_lst = [
        pdb_id_col,
        pdb_code_col,
        chainid_col,
        swiss_id_col,
        prot_col,
        method_col,
        resolution_col,
        r_factor_col,
        seq_col,
        len_col,
    ]

    record_dict = {}

    for col in col_lst:
        record_dict[col] = list()

    fix_lst = list()
    if fix_dict is not None:
        fix_lst = list(fix_dict.keys())

    for record in tqdm(
        list(pdbaa_dict.keys()), desc="Searching pdbaa", position=0, leave=True
    ):

        record = pdbaa_dict[record]

        desc = get_record_desc(record)

        swiss_id = get_pdbaa_swiss_id(desc)

        pdb_code = get_pdbaa_pdb_code(desc)
        chainid = get_pdbaa_chainid(desc)
        pdb_id = get_pdb_id(pdb_code, chainid)

        get_seq = False

        if search_col == swiss_id_col:
            if swiss_id in search_lst:
                get_seq = True

        elif search_col == pdb_id_col:
            if pdb_id in search_lst:
                get_seq = True

        if len(fix_lst) > 0:
            if not get_seq:
                if pdb_id in fix_lst:
                    fix_lst.remove(pdb_id)
                    swiss_id = fix_dict[pdb_id]
                    get_seq = True

        if get_seq:

            seq = get_record_seq(record)
            length = len(seq)

            get_lst = True
            if min_length is not None:
                if length <= min_length:
                    get_lst = False

            if get_lst:

                resolution = get_pdbaa_resolution(desc)

                add_dict = True
                if min_resolution is not None:
                    if resolution <= min_resolution:
                        add_dict = False

                if add_dict:

                    record_dict[pdb_id_col].append(pdb_id)
                    record_dict[pdb_code_col].append(pdb_code)
                    record_dict[chainid_col].append(chainid)
                    record_dict[method_col].append(get_pdbaa_method(desc))
                    record_dict[resolution_col].append(resolution)
                    record_dict[r_factor_col].append(get_pdbaa_r_factor(desc))
                    record_dict[swiss_id_col].append(swiss_id)
                    record_dict[prot_col].append(get_pdbaa_prot(desc))
                    record_dict[seq_col].append(seq)
                    record_dict[len_col].append(length)

    df = pd.DataFrame(record_dict)

    total_entry = len(lst_col(df, pdb_code_col, unique=True))
    total_chains = len(lst_col(df, pdb_id_col, unique=True))

    print(
        f"Your search identified {total_chains} chains from {total_entry} PDB entries!"
    )

    if entry_table_path is not None:
        save_table(entry_table_path, df)
    else:
        return df
