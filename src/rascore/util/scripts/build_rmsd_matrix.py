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
import itertools
import concurrent.futures

from ..functions.coord import load_coord, get_seq_lst, sup_without_map, calc_rmsd
from ..functions.table import fix_val, order_rows, lst_col, merge_dicts
from ..functions.lst import res_to_lst, lst_unique
from ..functions.path import save_matrix, save_json
from ..functions.col import core_path_col, modelid_col, chainid_col


def add_index_dict(final_dict, index_dict):

    for coord_path in list(index_dict.keys()):
        if coord_path not in list(final_dict.keys()):
            final_dict[coord_path] = dict()
        for modelid in list(index_dict[coord_path].keys()):
            if modelid not in list(final_dict[coord_path].keys()):
                final_dict[coord_path][modelid] = dict()
            for chainid in list(index_dict[coord_path][modelid].keys()):
                final_dict[coord_path][modelid][chainid] = index_dict[coord_path][
                    modelid
                ][chainid]

    return final_dict


def build_coord_dict(coord_path):

    return {coord_path: load_coord(coord_path)}


def build_seq_dict(structure, coord_path, modelid, chainid):

    modelid = fix_val(modelid, return_int=True)

    return {coord_path: {modelid: {chainid: get_seq_lst(structure[modelid][chainid])}}}


def build_sup_dict(
    coord_path,
    modelid,
    chainid,
    ref_chain,
    mob_chain,
    ref_seq_lst=None,
    mob_seq_lst=None,
    sup_resid_lst=None,
    pair_aln=True,
):

    sup = sup_without_map(
        ref_chain,
        mob_chain,
        ref_seq_lst=ref_seq_lst,
        mob_seq_lst=mob_seq_lst,
        sup_resids=sup_resid_lst,
        pair_aln=pair_aln,
    )

    return {coord_path: {str(modelid): {chainid: sup[1]}}}


def check_rmsd_dict(
    rmsd_dict, coord_path_1, modelid_1, chainid_1, coord_path_2, modelid_2, chainid_2
):

    rmsd_dist = None
    if coord_path_1 in list(rmsd_dict.keys()):
        if str(modelid_1) in list(rmsd_dict[coord_path_1].keys()):
            if chainid_1 in list(rmsd_dict[coord_path_1][str(modelid_1)].keys()):
                if coord_path_2 in list(
                    rmsd_dict[coord_path_1][str(modelid_1)][chainid_1].keys()
                ):
                    if modelid_2 in list(
                        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][
                            coord_path_2
                        ].keys()
                    ):
                        if chainid_2 in list(
                            rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][
                                coord_path_2
                            ][str(modelid_2)].keys()
                        ):
                            rmsd_dist = rmsd_dict[coord_path_1][str(modelid_1)][
                                chainid_1
                            ][coord_path_2][str(modelid_2)][chainid_2]

    return rmsd_dist


def append_rmsd_dict(
    rmsd_dict,
    rmsd_dist,
    coord_path_1,
    modelid_1,
    chainid_1,
    coord_path_2,
    modelid_2,
    chainid_2,
):

    if coord_path_1 not in list(rmsd_dict.keys()):
        rmsd_dict[coord_path_1] = dict()
    if str(modelid_1) not in list(rmsd_dict[coord_path_1].keys()):
        rmsd_dict[coord_path_1][str(modelid_1)] = dict()
    if chainid_1 not in list(rmsd_dict[coord_path_1][str(modelid_1)].keys()):
        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1] = dict()
    if coord_path_2 not in list(
        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1].keys()
    ):
        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][coord_path_2] = dict()
    if str(modelid_2) not in list(
        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][coord_path_2].keys()
    ):
        rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][coord_path_2][
            str(modelid_2)
        ] = dict()

    rmsd_dict[coord_path_1][str(modelid_1)][chainid_1][coord_path_2][str(modelid_2)][
        chainid_2
    ] = rmsd_dist

    return rmsd_dict


def prep_rmsd_job(
    i_df,
    j_df,
    i_index,
    j_index,
    coord_dict,
    seq_dict=None,
    sup_resid_lst=None,
    rmsd_resid_lst=None,
    rmsd_dict=None,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    i_coord_path = i_df.at[i_index, coord_path_col]
    i_modelid = i_df.at[i_index, modelid_col]
    i_chainid = i_df.at[i_index, chainid_col]

    j_coord_path = j_df.at[j_index, coord_path_col]
    j_modelid = j_df.at[j_index, modelid_col]
    j_chainid = j_df.at[j_index, chainid_col]

    i_chain = coord_dict[i_coord_path][fix_val(i_modelid, return_int=True)][i_chainid]
    j_chain = coord_dict[j_coord_path][fix_val(j_modelid, return_int=True)][j_chainid]

    if seq_dict is not None:
        i_seq_lst = seq_dict[i_coord_path][fix_val(i_modelid, return_int=True)][
            i_chainid
        ]
        j_seq_lst = seq_dict[j_coord_path][fix_val(j_modelid, return_int=True)][
            j_chainid
        ]
    else:
        i_seq_lst = None
        j_seq_lst = None

    rmsd_dist = None
    if rmsd_dict is not None:
        rmsd_dist = check_rmsd_dict(
            rmsd_dict,
            i_coord_path,
            i_modelid,
            i_chainid,
            j_coord_path,
            j_modelid,
            j_chainid,
        )

        if rmsd_dist is None:
            rmsd_dist = check_rmsd_dict(
                rmsd_dict,
                j_coord_path,
                j_modelid,
                j_chainid,
                i_coord_path,
                i_modelid,
                i_chainid,
            )

    return (
        i_chain,
        j_chain,
        i_seq_lst,
        j_seq_lst,
        sup_resid_lst,
        rmsd_resid_lst,
        rmsd_dist,
    )


def calc_rmsd_dist(
    i,
    j,
    i_chain,
    j_chain,
    i_seq_lst=None,
    j_seq_lst=None,
    sup_resids=None,
    rmsd_resids=None,
    rmsd_atomids="CA",
    pair_aln=True,
    rmsd_dist=None,
):

    if rmsd_dist is None:
        sup = sup_without_map(
            i_chain,
            j_chain,
            ref_seq_lst=i_seq_lst,
            mob_seq_lst=j_seq_lst,
            sup_resids=sup_resids,
            pair_aln=pair_aln,
        )

        rmsd_dist = calc_rmsd(
            sup[0],
            sup[1],
            map_dict=sup[2],
            rmsd_resids=rmsd_resids,
            rmsd_atomids=rmsd_atomids,
        )

    return (i, j, rmsd_dist)


def build_rmsd_matrix(
    fit_df,
    rmsd_matrix_path,
    pred_df=None,
    sup_resids=None,
    rmsd_resids=None,
    rmsd_atomids="CA",
    pair_aln=True,
    rmsd_dict=None,
    rmsd_json_path=None,
    coord_path_col=None,
    num_cpu=1,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    sup_resid_lst = None
    if sup_resids is not None:
        sup_resid_lst = res_to_lst(sup_resids)

    rmsd_resid_lst = None
    if rmsd_resids is not None:
        rmsd_resid_lst = res_to_lst(rmsd_resids)

    fit_df = order_rows(fit_df)

    j_df = fit_df.copy(deep=True)

    coord_path_lst = lst_col(fit_df, coord_path_col, unique=True)

    coord_df = fit_df.copy(deep=True)

    if pred_df is None:
        i_df = fit_df.copy(deep=True)
    else:
        pred_df = order_rows(pred_df)
        i_df = pred_df.copy(deep=True)

        coord_path_lst += lst_col(pred_df, coord_path_col, unique=True)
        coord_df = pd.concat([coord_df, pred_df])

    i_index_lst = list(i_df.index.values)
    j_index_lst = list(j_df.index.values)

    if pred_df is not None:
        index_pairs = itertools.product(i_index_lst, j_index_lst)
    else:
        index_pairs = itertools.combinations(i_index_lst, 2)

    coord_df = coord_df.reset_index(drop=True)

    coord_path_lst = lst_unique(coord_path_lst)

    coord_dict = dict()

    if num_cpu == 1:

        for coord_path in tqdm(
            coord_path_lst, desc="Loading coordinates", position=0, leave=True
        ):
            coord_dict = merge_dicts([coord_dict, build_coord_dict(coord_path)])

    else:

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:

            job_lst = [
                executor.submit(build_coord_dict, coord_path)
                for coord_path in coord_path_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Loading coordinates",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                coord_dict = merge_dicts([coord_dict, job.result()])

    seq_dict = None

    if pair_aln:

        seq_dict = dict()

        if num_cpu == 1:

            for index in tqdm(
                list(coord_df.index.values),
                desc="Loading sequences",
                position=0,
                leave=True,
            ):

                coord_path = coord_df.at[index, coord_path_col]
                modelid = coord_df.at[index, modelid_col]
                chainid = coord_df.at[index, chainid_col]

                structure = coord_dict[coord_path]

                seq_dict = add_index_dict(
                    seq_dict,
                    build_seq_dict(
                        structure,
                        coord_path,
                        modelid,
                        chainid,
                    ),
                )

        else:

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=num_cpu
            ) as executor:

                job_lst = list()

                for index in list(coord_df.index.values):

                    coord_path = coord_df.at[index, coord_path_col]
                    modelid = coord_df.at[index, modelid_col]
                    chainid = coord_df.at[index, chainid_col]

                    structure = coord_dict[coord_path]

                    job_lst.append(
                        executor.submit(
                            build_seq_dict,
                            structure,
                            coord_path,
                            modelid,
                            chainid,
                        )
                    )

                for job in tqdm(
                    concurrent.futures.as_completed(job_lst),
                    desc="Loading sequences",
                    total=len(job_lst),
                    miniters=1,
                    position=0,
                    leave=True,
                ):

                    seq_dict = add_index_dict(seq_dict, job.result())

    matrix = np.zeros((len(i_index_lst), len(j_index_lst)))

    if num_cpu == 1:

        for i_index, j_index in tqdm(
            index_pairs, desc="Building RMSD matrix", position=0, leave=True
        ):

            (
                i_chain,
                j_chain,
                i_seq_lst,
                j_seq_lst,
                sup_resid_lst,
                rmsd_resid_lst,
                rmsd_dist,
            ) = prep_rmsd_job(
                i_df,
                j_df,
                i_index,
                j_index,
                coord_dict,
                seq_dict=seq_dict,
                sup_resid_lst=sup_resid_lst,
                rmsd_resid_lst=rmsd_resid_lst,
                rmsd_dict=rmsd_dict,
                coord_path_col=coord_path_col,
            )

            result = calc_rmsd_dist(
                i_index,
                j_index,
                i_chain,
                j_chain,
                i_seq_lst=i_seq_lst,
                j_seq_lst=j_seq_lst,
                sup_resids=sup_resid_lst,
                rmsd_resids=rmsd_resid_lst,
                rmsd_atomids=rmsd_atomids,
                pair_aln=pair_aln,
                rmsd_dist=rmsd_dist,
            )

            i_index = result[0]
            j_index = result[1]
            dist = result[2]

            matrix[i_index, j_index] = dist

            if pred_df is None:
                matrix[j_index, i_index] = dist

    else:

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:

            job_lst = list()
            for i_index, j_index in index_pairs:

                (
                    i_chain,
                    j_chain,
                    i_seq_lst,
                    j_seq_lst,
                    sup_resid_lst,
                    rmsd_resid_lst,
                    rmsd_dist,
                ) = prep_rmsd_job(
                    i_df,
                    j_df,
                    i_index,
                    j_index,
                    coord_dict,
                    seq_dict=seq_dict,
                    sup_resid_lst=sup_resid_lst,
                    rmsd_resid_lst=rmsd_resid_lst,
                    rmsd_dict=rmsd_dict,
                    coord_path_col=coord_path_col,
                )

                job_lst.append(
                    executor.submit(
                        calc_rmsd_dist,
                        i_index,
                        j_index,
                        i_chain,
                        j_chain,
                        i_seq_lst=i_seq_lst,
                        j_seq_lst=j_seq_lst,
                        sup_resids=sup_resid_lst,
                        rmsd_resids=rmsd_resid_lst,
                        rmsd_atomids=rmsd_atomids,
                        pair_aln=pair_aln,
                        rmsd_dist=rmsd_dist,
                    )
                )

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Building RMSD matrix",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                result = job.result()

                i_index = result[0]
                j_index = result[1]
                dist = result[2]

                matrix[i_index, j_index] = dist

                if pred_df is None:
                    matrix[j_index, i_index] = dist

    save_matrix(rmsd_matrix_path, matrix)

    if rmsd_json_path is not None:

        if rmsd_dict is None:
            rmsd_dict = dict()

        if pred_df is not None:
            index_pairs = itertools.product(i_index_lst, j_index_lst)
        else:
            index_pairs = itertools.combinations(i_index_lst, 2)

        if num_cpu == 1:

            for i_index, j_index in tqdm(
                index_pairs, desc="Building RMSD dictionary", position=0, leave=True
            ):
                i_coord_path = i_df.at[i_index, coord_path_col]
                i_modelid = i_df.at[i_index, modelid_col]
                i_chainid = i_df.at[i_index, chainid_col]

                j_coord_path = j_df.at[j_index, coord_path_col]
                j_modelid = j_df.at[j_index, modelid_col]
                j_chainid = j_df.at[j_index, chainid_col]

                rmsd_dist = matrix[i_index, j_index]

                rmsd_dict = append_rmsd_dict(
                    rmsd_dict,
                    rmsd_dist,
                    i_coord_path,
                    i_modelid,
                    i_chainid,
                    j_coord_path,
                    j_modelid,
                    j_chainid,
                )

        else:

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=num_cpu
            ) as executor:

                if pred_df is not None:
                    index_pairs = itertools.product(i_index_lst, j_index_lst)
                else:
                    index_pairs = itertools.combinations(i_index_lst, 2)

                job_lst = list()

                for i_index, j_index in index_pairs:
                    i_coord_path = i_df.at[i_index, coord_path_col]
                    i_modelid = i_df.at[i_index, modelid_col]
                    i_chainid = i_df.at[i_index, chainid_col]

                    j_coord_path = j_df.at[j_index, coord_path_col]
                    j_modelid = j_df.at[j_index, modelid_col]
                    j_chainid = j_df.at[j_index, chainid_col]

                    rmsd_dist = matrix[i_index, j_index]

                    job_lst.append(
                        executor.submit(
                            append_rmsd_dict,
                            rmsd_dict,
                            rmsd_dist,
                            i_coord_path,
                            i_modelid,
                            i_chainid,
                            j_coord_path,
                            j_modelid,
                            j_chainid,
                        )
                    )

                for job in tqdm(
                    concurrent.futures.as_completed(job_lst),
                    desc="Building RMSD dictionary",
                    total=len(job_lst),
                    miniters=1,
                    position=0,
                    leave=True,
                ):

                    rmsd_dict = job.result()

        save_json(rmsd_json_path, rmsd_dict)

    print("Built RMSD matrix!")
