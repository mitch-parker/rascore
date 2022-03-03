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

from tqdm import tqdm
import concurrent.futures

from ..functions.path import save_json
from ..functions.table import merge_dicts
from ..functions.dih import calc_bb_angle, calc_sc_angle
from ..functions.coord import (
    load_coord,
    get_modelid,
    get_chainid,
    get_resid,
    get_resname,
    get_resid_str,
    is_aa,
)
from ..functions.col import resname_col, bb_col_lst, sc_col_lst
from ..functions.lst import type_lst
from ..functions.table import merge_dicts
from ..functions.path import save_json


def build_dih_dict(coord_path):

    structure = load_coord(coord_path)

    dih_dict = {coord_path: dict()}

    for model in structure:
        modelid = get_modelid(model)
        dih_dict[coord_path][modelid] = dict()
        for chain in model:
            chainid = get_chainid(chain)

            dih_dict[coord_path][modelid][chainid] = dict()

            resid_lst = [get_resid(residue) for residue in chain if is_aa(residue)]
            resname_lst = [get_resname(residue) for residue in chain if is_aa(residue)]
            resstr_lst = [get_resid_str(residue) for residue in chain if is_aa(residue)]

            for index, resid in enumerate(resid_lst):

                curr_resid = resid
                curr_resname = resname_lst[index]
                curr_resid_str = resstr_lst[index]

                dih_dict[coord_path][modelid][chainid][curr_resid_str] = dict()
                dih_dict[coord_path][modelid][chainid][curr_resid_str][
                    resname_col
                ] = curr_resname

                if index != 0 and index != len(resid_lst) - 1:
                    prev_resid = resid_lst[index - 1]
                    next_resid = resid_lst[index + 1]

                    for bb_col in bb_col_lst:

                        bb_val = calc_bb_angle(
                            structure=structure,
                            modelid=modelid,
                            chainid=chainid,
                            curr_resid=curr_resid,
                            prev_resid=prev_resid,
                            next_resid=next_resid,
                            angle=bb_col,
                        )

                        dih_dict[coord_path][modelid][chainid][curr_resid_str][
                            bb_col
                        ] = bb_val

                else:
                    for bb_col in bb_col_lst:
                        dih_dict[coord_path][modelid][chainid][curr_resid_str][
                            bb_col
                        ] = 999.00

                for sc_col in sc_col_lst:

                    sc_val = calc_sc_angle(
                        structure=structure,
                        modelid=modelid,
                        chainid=chainid,
                        resid=curr_resid,
                        angle=sc_col,
                    )
                    dih_dict[coord_path][modelid][chainid][curr_resid_str][
                        sc_col
                    ] = sc_val

    return dih_dict


def prep_dih(coord_paths, dih_json_path=None, num_cpu=1):

    coord_path_lst = type_lst(coord_paths)

    dih_dict = dict()

    if num_cpu == 1:
        for coord_path in tqdm(
            coord_path_lst,
            desc="Preparing dihedrals",
            position=0,
            leave=True,
        ):
            dih_dict = merge_dicts([dih_dict, build_dih_dict(coord_path)])

    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(build_dih_dict, coord_path)
                for coord_path in coord_path_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Preparing dihedrals",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):
                dih_dict = merge_dicts([dih_dict, job.result()])

    print("Prepared dihedrals!")

    if dih_json_path is not None:
        save_json(dih_json_path, dih_dict)
    else:
        return dih_dict
