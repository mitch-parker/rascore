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
import concurrent.futures

from ..functions.seq import calc_seq_id, load_record_lst, get_record_seq
from ..functions.lst import lst_to_str, res_to_lst, lst_nums, str_to_lst
from ..functions.coord import load_coord, get_resname, has_resid
from ..functions.col import (
    core_path_col,
    modelid_col,
    chainid_col,
    seq_col,
    mut_status_col,
    mut_pos_col,
    uniprot_id_col,
)
from ..functions.table import get_df_at_index, fix_val
from ..functions.path import save_table, get_seq_path
from ..functions.url import uniprot_url
from ..functions.download import download_file


def build_mut_df(df, index, uniprot_dict, resid_lst=None, coord_path_col=None):

    if coord_path_col is None:
        coord_path_col = core_path_col

    df = get_df_at_index(df, index)

    coord_path = df.at[index, coord_path_col]
    modelid = df.at[index, modelid_col]
    chainid = df.at[index, chainid_col]

    structure = load_coord(coord_path)

    uniprot_acc_lst = list(uniprot_dict.keys())

    mut_status_dict = dict()
    mut_pos_dict = dict()
    uniprot_id_lst = list()

    for uniprot_acc in uniprot_acc_lst:

        mut_status_dict[uniprot_acc] = list()
        mut_pos_dict[uniprot_acc] = list()

        uniprot_seq = uniprot_dict[uniprot_acc][seq_col]

        coord_seq = ""
        for resid in lst_nums(1, len(uniprot_seq)):
            resname = "-"
            if has_resid(structure, chainid, resid, modelid=modelid):
                if resname != "X":
                    resname = get_resname(
                        structure[fix_val(modelid, return_int=True)][chainid][resid],
                        letter=True,
                    )
                    ref = uniprot_dict[uniprot_acc][resid]

                    if resname != ref:
                        if resid is not None:
                            add_mut = True
                            if resid_lst is not None:
                                if resid not in resid_lst:
                                    add_mut = False
                            if add_mut:
                                mut_pos = f"{ref}{resid}"
                                mut_status_dict[uniprot_acc].append(mut_pos + resname)
                                mut_pos_dict[uniprot_acc].append(mut_pos)

            coord_seq += resname

        uniprot_id_lst.append(calc_seq_id(coord_seq, uniprot_seq, aln=False))

    uniprot_acc = uniprot_acc_lst[uniprot_id_lst.index(max(uniprot_id_lst))]

    mut_status_lst = mut_status_dict[uniprot_acc]
    mut_pos_lst = mut_pos_dict[uniprot_acc]

    df.at[index, mut_status_col] = lst_to_str(mut_status_lst, empty="WT")
    df.at[index, mut_pos_col] = lst_to_str(mut_pos_lst, empty="WT")
    df.at[index, uniprot_id_col] = uniprot_acc

    return df


def annot_mut(
    df,
    uniprot_accs,
    mut_table_path=None,
    resids=None,
    seq_dir=None,
    coord_path_col=None,
    num_cpu=1,
):

    df = df.reset_index(drop=True)

    if type(uniprot_accs) == list:
        uniprot_acc_lst = uniprot_accs
    else:
        uniprot_acc_lst = str_to_lst(uniprot_accs, sep_txt=" ")

    if resids is not None:
        resid_lst = res_to_lst(resids)
    else:
        resid_lst = None

    uniprot_dict = dict()

    for uniprot_acc in uniprot_acc_lst:

        uniprot_dict[uniprot_acc] = dict()

        fasta_url = f"{uniprot_url}{uniprot_acc}.fasta"

        fasta_file = get_seq_path(uniprot_acc, dir_path=seq_dir)

        download_file(fasta_url, fasta_file)

        for record in load_record_lst(fasta_file):

            seq = get_record_seq(record)

        uniprot_dict[uniprot_acc][seq_col] = seq

        for i, resname in enumerate(seq):

            uniprot_dict[uniprot_acc][i + 1] = resname

    mut_df = pd.DataFrame()

    if num_cpu == 1:
        for index in tqdm(
            list(df.index.values), desc="Annotating mutations", position=0, leave=True
        ):
            mut_df = pd.concat(
                [
                    mut_df,
                    build_mut_df(
                        df,
                        index,
                        uniprot_dict,
                        resid_lst=resid_lst,
                        coord_path_col=coord_path_col,
                    ),
                ],
                sort=False,
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    build_mut_df,
                    df,
                    index,
                    uniprot_dict,
                    resid_lst=resid_lst,
                )
                for index in list(df.index.values)
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Annotating mutations",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                mut_df = pd.concat([mut_df, job.result()], sort=False)

    mut_df = mut_df.reset_index(drop=True)

    print("Annotated mutations!")

    if mut_table_path is not None:
        save_table(mut_table_path, mut_df)
    else:
        return mut_df
