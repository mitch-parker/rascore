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

from ..functions.coord import load_coord, calc_atom_dist, wmhb_name, no_hb_name
from ..functions.lst import type_lst, lst_nums, lst_to_str, str_to_lst
from ..functions.col import (
    get_dist_col,
    core_path_col,
    modelid_col,
    chainid_col,
    vect_1_col,
    vect_2_col,
    hb_status_col,
    hb_angle_1_col,
    hb_angle_2_col,
    wmhb_angle_col,
    outlier_col,
)
from ..functions.table import get_df_at_index
from ..functions.path import modify_coord_path, save_table


def get_index_dist(
    df,
    index,
    x_resids,
    y_resids,
    x_atomids=None,
    y_atomids=None,
    shared_resids=None,
    shared_atomids=None,
    atom_dist_col_lst=None,
    vect_1_col_lst=None,
    vect_2_col_lst=None,
    hb_status_col_lst=None,
    hb_angle_1_col_lst=None,
    hb_angle_2_col_lst=None,
    wmhb_angle_col_lst=None,
    outlier_col_lst=None,
    check_hb=False,
    use_h=False,
    hb_sc=True,
    hb_bb=True,
    min_hb_dist=2.0,
    max_hb_dist=3.2,
    min_wmhb_dist=2.0,
    max_wmhb_dist=3.0,
    min_hb_angle=90,
    max_hb_angle=180,
    min_wmhb_angle=80,
    max_wmhb_angle=140,
    coord_path_col=None,
):

    if coord_path_col is None:
        coord_path_col = core_path_col

    index_df = get_df_at_index(df, index)

    coord_path = index_df.at[index, coord_path_col]
    modelid = index_df.at[index, modelid_col]
    chainid = index_df.at[index, chainid_col]

    val_lst = [x_resids, y_resids]

    if x_atomids is not None:
        val_lst.append(x_atomids)
    if y_atomids is not None:
        val_lst.append(y_atomids)

    if shared_resids is not None:
        val_lst.append(shared_resids)
    if shared_atomids is not None:
        val_lst.append(shared_atomids)

    for i, val in enumerate(val_lst):
        val_lst[i] = type_lst(val)

    x_lst = lst_nums(0, len(val_lst[0]) - 1)

    if use_h:
        coord_path = modify_coord_path(coord_path, return_pdb=True, add_h=True)

    structure = load_coord(coord_path)

    for x in x_lst:

        resid_1 = val_lst[0][x]
        resid_2 = val_lst[1][x]

        if x_atomids is not None:
            atomid_1 = val_lst[2][x]
        else:
            atomid_1 = None

        if y_atomids is not None:
            atomid_2 = val_lst[3][x]
        else:
            atomid_2 = None

        if atom_dist_col_lst is None:
            cont_atom_dist_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
            )
        else:
            cont_atom_dist_col = atom_dist_col_lst[x]

        if vect_1_col_lst is None:
            cont_vect_1_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=vect_1_col,
            )
        else:
            cont_vect_1_col = vect_1_col_lst[x]

        if vect_2_col_lst is None:
            cont_vect_2_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=vect_2_col,
            )
        else:
            cont_vect_2_col = vect_2_col_lst[x]

        if hb_status_col_lst is None:
            cont_hb_status_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=hb_status_col,
            )
        else:
            cont_hb_status_col = hb_status_col_lst[x]

        if hb_angle_1_col_lst is None:
            cont_hb_angle_1_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=hb_angle_1_col,
            )
        else:
            cont_hb_angle_1_col = hb_angle_1_col_lst[x]

        if hb_angle_2_col_lst is None:
            cont_hb_angle_2_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=hb_angle_2_col,
            )
        else:
            cont_hb_angle_2_col = hb_angle_2_col_lst[x]

        if wmhb_angle_col_lst is None:
            cont_wmhb_angle_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=wmhb_angle_col,
            )
        else:
            cont_wmhb_angle_col = wmhb_angle_col_lst[x]

        if outlier_col_lst is None:
            cont_outlier_col = get_dist_col(
                resid_1,
                resid_2,
                x_atomid=lst_to_str(atomid_1, join_txt="+"),
                y_atomid=lst_to_str(atomid_2, join_txt="+"),
                ext=outlier_col,
            )
        else:
            cont_outlier_col = outlier_col_lst[x]

        resid_lst = [resid_1, resid_2]

        for i, resid in enumerate(resid_lst):
            if type(resid) == str:
                resid = str_to_lst(df.at[index, resid])[0]
            resid_lst[i] = resid

        result = calc_atom_dist(
            structure=structure,
            chainid_1=chainid,
            modelid_1=modelid,
            resid_1=resid_lst[0],
            atomid_1=atomid_1,
            chainid_2=chainid,
            modelid_2=modelid,
            resid_2=resid_lst[1],
            atomid_2=atomid_2,
            check_hb=check_hb,
            use_h=use_h,
            hb_sc=hb_sc,
            hb_bb=hb_bb,
            min_hb_dist=min_hb_dist,
            max_hb_dist=max_hb_dist,
            min_wmhb_dist=min_wmhb_dist,
            max_wmhb_dist=max_wmhb_dist,
            min_hb_angle=min_hb_angle,
            max_hb_angle=max_hb_angle,
            min_wmhb_angle=min_wmhb_angle,
            max_wmhb_angle=max_wmhb_angle,
            return_vect=True,
        )

        atom_dist = result[0]
        vect_1 = result[1]
        vect_2 = result[2]

        if check_hb:
            hb_status = result[3]
            hb_angle_1 = result[4]
            hb_angle_2 = result[5]
            wmhb_angle = result[6]
            outlier_status = result[7]

            if hb_status == wmhb_name:
                if shared_resids is not None:
                    shared_resid = val_lst[4][x]

                    if shared_atomids is not None:
                        shared_atomid = val_lst[5][x]
                    else:
                        shared_atomid = None

                    shared_result_1 = calc_atom_dist(
                        structure=structure,
                        chainid_1=chainid,
                        modelid_1=modelid,
                        resid_1=resid_lst[0],
                        atomid_1=atomid_1,
                        chainid_2=chainid,
                        modelid_2=modelid,
                        resid_2=shared_resid,
                        atomid_2=shared_atomid,
                        check_hb=check_hb,
                        use_h=use_h,
                        hb_sc=hb_sc,
                        hb_bb=hb_bb,
                        min_hb_dist=min_hb_dist,
                        max_hb_dist=max_hb_dist,
                        min_wmhb_dist=min_wmhb_dist,
                        max_wmhb_dist=max_wmhb_dist,
                        min_hb_angle=min_hb_angle,
                        max_hb_angle=max_hb_angle,
                        min_wmhb_angle=min_wmhb_angle,
                        max_wmhb_angle=max_wmhb_angle,
                        return_vect=True,
                    )
                    shared_result_2 = calc_atom_dist(
                        structure=structure,
                        chainid_1=chainid,
                        modelid_1=modelid,
                        resid_1=shared_resid,
                        atomid_1=shared_atomid,
                        chainid_2=chainid,
                        modelid_2=modelid,
                        resid_2=resid_lst[1],
                        atomid_2=atomid_2,
                        check_hb=check_hb,
                        use_h=use_h,
                        hb_sc=hb_sc,
                        hb_bb=hb_bb,
                        min_hb_dist=min_hb_dist,
                        max_hb_dist=max_hb_dist,
                        min_wmhb_dist=min_wmhb_dist,
                        max_wmhb_dist=max_wmhb_dist,
                        min_hb_angle=min_hb_angle,
                        max_hb_angle=max_hb_angle,
                        min_wmhb_angle=min_wmhb_angle,
                        max_wmhb_angle=max_wmhb_angle,
                        return_vect=True,
                    )

                    if (shared_result_1[3] == no_hb_name) or (
                        shared_result_2[3] == no_hb_name
                    ):
                        hb_status = no_hb_name
                        outlier_status = True

        index_df.at[
            index,
            cont_atom_dist_col,
        ] = atom_dist

        index_df.at[
            index,
            cont_vect_1_col,
        ] = vect_1

        index_df.at[
            index,
            cont_vect_2_col,
        ] = vect_2

        if check_hb:
            index_df.at[
                index,
                cont_hb_status_col,
            ] = hb_status

            index_df.at[
                index,
                cont_hb_angle_1_col,
            ] = hb_angle_1

            index_df.at[
                index,
                cont_hb_angle_2_col,
            ] = hb_angle_2

            index_df.at[
                index,
                cont_wmhb_angle_col,
            ] = wmhb_angle
            index_df.at[
                index,
                cont_outlier_col,
            ] = outlier_status

    return index_df


def build_dist_table(
    df,
    x_resids,
    y_resids,
    x_atomids=None,
    y_atomids=None,
    shared_resids=None,
    shared_atomids=None,
    dist_table_path=None,
    atom_dist_col_lst=None,
    hb_status_col_lst=None,
    hb_angle_1_col_lst=None,
    hb_angle_2_col_lst=None,
    wmhb_angle_col_lst=None,
    outlier_col_lst=None,
    check_hb=False,
    use_h=False,
    hb_sc=True,
    hb_bb=True,
    min_hb_dist=2.0,
    max_hb_dist=3.2,
    min_wmhb_dist=2.0,
    max_wmhb_dist=3.0,
    min_hb_angle=90,
    max_hb_angle=180,
    min_wmhb_angle=80,
    max_wmhb_angle=140,
    coord_path_col=None,
):

    df = df.reset_index(drop=True)

    dist_df = pd.DataFrame()

    for index in tqdm(
        list(df.index.values), desc="Building distance table", position=0, leave=True
    ):
        dist_df = pd.concat(
            [
                dist_df,
                get_index_dist(
                    df,
                    index,
                    x_resids,
                    y_resids,
                    x_atomids=x_atomids,
                    y_atomids=y_atomids,
                    shared_resids=shared_resids,
                    shared_atomids=shared_atomids,
                    atom_dist_col_lst=atom_dist_col_lst,
                    hb_status_col_lst=hb_status_col_lst,
                    hb_angle_1_col_lst=hb_angle_1_col_lst,
                    hb_angle_2_col_lst=hb_angle_2_col_lst,
                    wmhb_angle_col_lst=wmhb_angle_col_lst,
                    outlier_col_lst=outlier_col_lst,
                    check_hb=check_hb,
                    use_h=use_h,
                    hb_sc=hb_sc,
                    hb_bb=hb_bb,
                    min_hb_dist=min_hb_dist,
                    max_hb_dist=max_hb_dist,
                    min_wmhb_dist=min_wmhb_dist,
                    max_wmhb_dist=max_wmhb_dist,
                    min_hb_angle=min_hb_angle,
                    max_hb_angle=max_hb_angle,
                    min_wmhb_angle=min_wmhb_angle,
                    max_wmhb_angle=max_wmhb_angle,
                    coord_path_col=coord_path_col,
                ),
            ],
            sort=False,
        )

    dist_df = dist_df.reset_index(drop=True)

    print("Built distance table!")

    if dist_table_path is not None:
        save_table(dist_table_path, dist_df)
    else:
        return dist_df
