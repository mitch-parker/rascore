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

from ..functions.lst import str_to_lst, lst_to_str, sort_lst, res_to_lst
from ..functions.table import extract_int, fix_val, get_df_at_index, fix_col
from ..functions.coord import (
    load_coord,
    get_neighbors,
    has_resid,
    resid_to_tuple,
    get_resname,
    get_residue_cont,
    get_reschainid,
    get_resnum,
    is_aa,
    calc_atom_dist,
)
from ..functions.col import (
    interf_path_col,
    atomid_cont_col,
    renum_path_col,
    modelid_col,
    chainid_col,
    interf_area_col,
    iso_col,
    interf_col,
    interf_cont_col,
    cb_dist_col,
    interf_id_col,
    pdb_id_col,
)
from ..functions.interf import calc_q_score
from ..functions.path import save_table


def calc_interf_cont_dist(
    coord_path,
    chainid,
    interf,
    interf_dict,
    cont_dist=5,
    modelid=None,
    interf_resid_lst=None,
):

    if modelid is None:
        modelid = 0

    interf_path = interf_dict[coord_path][chainid][interf][interf_path_col]

    structure = load_coord(interf_path)

    neighbors = get_neighbors(structure)

    resid_lst = str_to_lst(interf_dict[coord_path][chainid][interf][atomid_cont_col])

    interf_cont_lst = list()
    cb_dist_lst = list()

    add_interf = True
    if interf_resid_lst is not None:
        add_interf = False

    for resid in resid_lst:

        if interf_resid_lst is not None:
            if extract_int(resid) in interf_resid_lst:
                add_interf = True

        if resid is not None and add_interf:

            if has_resid(structure, chainid, resid_to_tuple(resid), modelid=modelid):
                residue = structure[fix_val(modelid, return_int=True)][chainid][
                    resid_to_tuple(resid)
                ]

                resname = get_resname(residue)

                if resname == "GLY":
                    atomid = "CA"
                else:
                    atomid = "CB"

                for cont_residue in get_residue_cont(
                    neighbors, residue, max_dist=cont_dist, level="R"
                ):
                    cont_chainid = get_reschainid(cont_residue)
                    cont_resid = get_resnum(cont_residue)
                    cont_resname = get_resname(cont_residue)

                    if is_aa(cont_residue) and cont_chainid != chainid:

                        if cont_resname == "GLY":
                            cont_atomid = "CA"
                        else:
                            cont_atomid = "CB"

                        if cont_resid is not None:

                            interf_cont = lst_to_str(
                                sort_lst([resid, cont_resid]),
                                join_txt=":",
                            )

                            if interf_cont not in interf_cont_lst:

                                cb_dist = calc_atom_dist(
                                    structure,
                                    modelid_1=modelid,
                                    chainid_1=chainid,
                                    resid_1=resid,
                                    atomid_1=atomid,
                                    modelid_2=modelid,
                                    chainid_2=cont_chainid,
                                    resid_2=cont_resid,
                                    atomid_2=cont_atomid,
                                )

                                interf_cont_lst.append(interf_cont)
                                cb_dist_lst.append(cb_dist)

    return interf_cont_lst, cb_dist_lst


def get_index_interf(
    df,
    index,
    interf_dict,
    interf_resid_lst=None,
    min_area=200,
    cont_dist=5,
    iso_interf=False,
    het_interf=False,
    search_interf_cont_lst=None,
    search_cb_dist_lst=None,
    search_max_dist=0.7,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = renum_path_col

    index_df = get_df_at_index(df, index)

    coord_path = index_df.at[index, coord_path_col]
    modelid = index_df.at[index, modelid_col]
    chainid = index_df.at[index, chainid_col]

    interf_df = pd.DataFrame()

    for interf in list(interf_dict[coord_path][chainid].keys()):

        add_interf = True

        if float(interf_dict[coord_path][chainid][interf][interf_area_col]) < min_area:
            add_interf = False

        interf_iso_status = interf_dict[coord_path][chainid][interf][iso_col]
        if iso_interf:
            if not interf_iso_status:
                add_interf = False
        if het_interf:
            if interf_iso_status:
                add_interf = False

        if add_interf:

            interf_cont_lst, cb_dist_lst = calc_interf_cont_dist(
                coord_path,
                chainid,
                interf,
                interf_dict,
                cont_dist=cont_dist,
                modelid=modelid,
                interf_resid_lst=interf_resid_lst,
            )

            if search_interf_cont_lst is not None and search_cb_dist_lst is not None:
                if (
                    calc_q_score(
                        i_cont_lst=search_interf_cont_lst,
                        j_cont_lst=interf_cont_lst,
                        i_dist_lst=search_cb_dist_lst,
                        j_dist_lst=cb_dist_lst,
                    )
                    > search_max_dist
                ):
                    add_interf = False

        if add_interf:
            temp_df = index_df.copy(deep=True)

            for col in list(interf_dict[coord_path][chainid][interf].keys()):
                temp_df.at[index, col] = interf_dict[coord_path][chainid][interf][col]

            temp_df.at[index, interf_col] = interf
            temp_df.at[index, interf_cont_col] = lst_to_str(interf_cont_lst)
            temp_df.at[index, cb_dist_col] = lst_to_str(cb_dist_lst)

            interf_df = pd.concat([interf_df, temp_df], sort=False)

    return interf_df


def build_interf_table(
    df,
    interf_dict,
    interf_table_path=None,
    interf_resids=None,
    min_area=200,
    cont_dist=5,
    iso_interf=False,
    het_interf=False,
    search_coord_path=None,
    search_modelid=None,
    search_chainid=None,
    search_interf=None,
    search_max_dist=0.7,
    coord_path_col=None,
):

    df = df.reset_index(drop=True)

    if coord_path_col is None:
        coord_path_col = renum_path_col

    interf_df = pd.DataFrame()

    interf_resid_lst = res_to_lst(interf_resids)

    search_interf_cont_lst = None
    search_cb_dist_lst = None

    if (
        search_coord_path is not None
        and search_chainid is not None
        and search_interf is not None
    ):

        try:
            search_interf_cont_lst, search_cb_dist_lst = calc_interf_cont_dist(
                search_coord_path,
                search_chainid,
                str(search_interf),
                interf_dict,
                cont_dist=cont_dist,
                modelid=search_modelid,
                interf_resid_lst=interf_resid_lst,
            )
        except:
            search_interf_cont_lst, search_cb_dist_lst = calc_interf_cont_dist(
                search_coord_path,
                search_chainid,
                int(search_interf),
                interf_dict,
                cont_dist=cont_dist,
                modelid=search_modelid,
                interf_resid_lst=interf_resid_lst,
            )

    for index in tqdm(
        list(df.index.values), desc="Building interface table", position=0, leave=True
    ):

        interf_df = pd.concat(
            [
                interf_df,
                get_index_interf(
                    df,
                    index,
                    interf_dict,
                    min_area=min_area,
                    cont_dist=cont_dist,
                    interf_resid_lst=interf_resid_lst,
                    iso_interf=iso_interf,
                    het_interf=het_interf,
                    search_interf_cont_lst=search_interf_cont_lst,
                    search_cb_dist_lst=search_cb_dist_lst,
                    search_max_dist=search_max_dist,
                    coord_path_col=coord_path_col,
                ),
            ],
            sort=False,
        )

    if len(interf_df) > 0:
        interf_df = interf_df.reset_index(drop=True)

        interf_df[interf_col] = interf_df[interf_col].map(str)
        interf_df = fix_col(interf_df, interf_col)

        interf_df[interf_id_col] = interf_df[pdb_id_col] + interf_df[interf_col]

    if interf_table_path is not None:
        save_table(interf_table_path, interf_df)
    else:
        return interf_df

    print("Built interface table!")
