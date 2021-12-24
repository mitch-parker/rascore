# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import pandas as pd
from tqdm import tqdm
import concurrent.futures

from util.data import get_df_at_index, fix_val
from util.lst import lst_to_str, lst_unique, res_to_lst, type_lst
from util.path import save_table
from util.coord import (
    load_coord,
    get_residues,
    is_het,
    get_resname,
    get_resnum,
    get_reschainid,
    is_aa,
    get_neighbors,
    get_residue_cont,
)

from util.col import (
    core_path_col,
    modelid_col,
    chainid_col,
    bio_lig_col,
    pharm_lig_col,
    pharm_lig_site_col,
    pharm_lig_match_col,
    bio_lig_cont_col,
    pharm_lig_cont_col,
)
from util.chem import is_lig_match, get_lig_simi

from util.lig import lig_class_dict, lig_col_lst


def build_lig_df(
    df,
    index,
    site_dict=None,
    match_dict=None,
    other_site="Other",
    other_match="Unclassified",
    none_site="None",
    none_match="None",
    max_site_dist=5,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    df = get_df_at_index(df, index)

    coord_path = df.at[index, coord_path_col]
    modelid = df.at[index, modelid_col]
    chainid = df.at[index, chainid_col]

    lig_lst_dict = dict()

    for lig_col in lig_col_lst:

        lig_lst_dict[lig_col] = list()

    chain = load_coord(coord_path)[fix_val(modelid, return_int=True)][chainid]
    neighbors = get_neighbors(chain)

    bio_cont_lst = list()
    pharm_cont_lst = list()
    site_lst = list()
    match_lst = list()

    for residue in get_residues(chain):
        if is_het(residue):
            resname = get_resname(residue)
            lig_class = False
            is_bio = False
            is_pharm = False
            for lig_col in lig_col_lst:
                if resname in lig_class_dict[lig_col]:
                    lig_class = True
                    if resname not in lig_lst_dict[lig_col]:
                        lig_lst_dict[lig_col].append(resname)
                        if lig_col == bio_lig_col:
                            is_bio = True
                        if lig_col == pharm_lig_col:
                            is_pharm = True

            if not lig_class:
                if resname not in lig_lst_dict[pharm_lig_col]:
                    lig_lst_dict[pharm_lig_col].append(resname)
                    is_pharm = True

            if is_bio or is_pharm:
                cont_lst = lst_unique(
                    [
                        get_resnum(x)
                        for x in get_residue_cont(
                            neighbors,
                            residue,
                            max_dist=max_site_dist,
                            level="R",
                        )
                        if (get_resnum(x) < 50000)
                        and (get_reschainid(x) == chainid)
                        and (is_aa(x))
                    ]
                )

            if is_bio:
                bio_cont_lst += cont_lst

            if is_pharm:
                pharm_cont_lst += cont_lst

                if site_dict is not None:

                    max_cont = 0
                    site_status = other_site
                    for site_name, site_cont_lst in site_dict.items():
                        site_cont = len(
                            [x for x in site_cont_lst if x in pharm_cont_lst]
                        )
                        if site_cont > max_cont:
                            site_status = site_name
                            max_cont = site_cont
                    site_lst.append(site_status)

                if match_dict is not None:
                    if site_status in list(match_dict.keys()):
                        match_status = other_match
                        match_name_lst = list()
                        match_simi_lst = list()
                        for match_name, query_lst in match_dict[site_status].items():
                            is_match = is_lig_match(
                                resname,
                                query_lst,
                            )
                            if is_match:
                                for query in query_lst:
                                    match_name_lst.append(match_name)
                                    match_simi_lst.append(get_lig_simi(resname, query))
                        if len(match_name_lst) > 0:
                            match_status = match_name_lst[
                                match_simi_lst.index(max(match_simi_lst))
                            ]
                        match_lst.append(match_status)

    for lig_col in lig_col_lst:
        df.at[index, lig_col] = lst_to_str(lig_lst_dict[lig_col])

    df.at[index, pharm_lig_site_col] = lst_to_str(lst_unique(site_lst), empty=none_site)
    df.at[index, pharm_lig_match_col] = lst_to_str(
        lst_unique(match_lst), empty=none_match
    )
    df.at[index, bio_lig_cont_col] = lst_to_str(lst_unique(bio_cont_lst))

    df.at[index, pharm_lig_cont_col] = lst_to_str(
        lst_unique(pharm_cont_lst), empty=none_match
    )

    return df


def annot_lig(
    df,
    lig_table_path=None,
    site_dict=None,
    match_dict=None,
    other_site="Other",
    other_match="Unclassified",
    none_site="None",
    none_match="None",
    max_site_dist=5,
    coord_path_col=None,
    num_cpu=1,
):

    if site_dict is not None:
        for site_name, site_cont in site_dict.items():
            site_dict[site_name] = res_to_lst(site_cont)

    if match_dict is not None:
        for site_name in list(match_dict.keys()):
            for match_name, query_lst in match_dict[site_name].items():
                match_dict[site_name][match_name] = type_lst(query_lst)

    lig_df = pd.DataFrame()

    if num_cpu == 1:
        for index in tqdm(
            list(df.index.values), desc="Annotating ligands", position=0, leave=True
        ):

            lig_df = pd.concat(
                [
                    lig_df,
                    build_lig_df(
                        df,
                        index,
                        site_dict=site_dict,
                        match_dict=match_dict,
                        other_site=other_site,
                        other_match=other_match,
                        none_site=none_site,
                        none_match=none_match,
                        max_site_dist=max_site_dist,
                        coord_path_col=coord_path_col,
                    ),
                ],
                sort=False,
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    build_lig_df,
                    df,
                    index,
                    site_dict=site_dict,
                    match_dict=match_dict,
                    other_site=other_site,
                    other_match=other_match,
                    none_site=none_site,
                    none_match=none_match,
                    max_site_dist=max_site_dist,
                    coord_path_col=coord_path_col,
                )
                for index in list(df.index.values)
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Annotating ligands",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                lig_df = pd.concat([lig_df, job.result()], sort=False)

    lig_df = lig_df.reset_index(drop=True)

    print("Annotated ligands!")

    if lig_table_path is not None:
        save_table(lig_table_path, lig_df)
    else:
        return lig_df
