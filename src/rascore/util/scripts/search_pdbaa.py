# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Mitchell Isaac Parker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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
