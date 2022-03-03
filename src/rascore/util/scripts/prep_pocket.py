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
import subprocess
import io

import pandas as pd

from tqdm import tqdm
import concurrent.futures
import pymol2

from ..functions.lst import (
    sort_lst,
    lst_unique,
    lst_to_str,
    type_lst,
)
from ..functions.table import (
    mask_equal,
    title_str,
    merge_dicts,
)
from ..functions.path import (
    delete_path,
    get_file_path,
    get_dir_path,
    get_file_name,
    copy_path,
    get_pocket_path,
    append_path,
    save_json,
    pocket_str,
)
from ..functions.coord import (
    load_coord,
    get_reschainid,
    get_residues,
    get_resnum,
    is_aa,
)
from ..functions.col import (
    pocket_path_col,
    pocket_cont_col,
    pocket_status_col,
    pocket_score_col,
    pocket_volume_col,
    pocket_type_col,
    pocket_lig_col,
    pharm_lig_col,
)
from ..functions.lig import lig_col_lst, lig_lst_dict

pocket_bound_name = "Bound"
pocket_unbound_name = "Unbound"

pocket_info_dict = {
    "Druggability Score": "drug_score",
    "Volume": "pock_vol",
    "Number of Alpha Spheres": "nb_AS",
}


def run_fpocket(
    coord_path,
    pocket_dir=None,
    chainid_lst=None,
):

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

    coord_name = get_file_name(coord_path)
    coord_name = coord_name.replace("_core", "")

    run_path = get_file_path(coord_name, dir_str=pocket_str, dir_path=pocket_dir)

    copy_path(coord_path, run_path)

    cmd_lst = [
        "fpocket",
        "-f",
        run_path,
        "-d",
    ]

    df = pd.read_csv(
        io.StringIO(
            subprocess.Popen(cmd_lst, stdout=subprocess.PIPE)
            .communicate()[0]
            .decode("utf-8")
        ),
        delim_whitespace=True,
        dtype=str,
    )

    if chainid_lst is not None:
        for chainid in chainid_lst:
            df = mask_equal(df, "name_chain_1", chainid)
            df = mask_equal(df, "name_chain_2", chainid)

    pocket_dict = {coord_path: dict()}

    obj = coord_name.replace(".pdb", "")

    run_dir_path = f"{pocket_dir_path}/pocket_{obj}_out"

    with pymol2.PyMOL() as pymol:

        cmd = pymol.cmd

        cmd.load(coord_path, obj)

        for i, index in enumerate(list(df.index.values)):

            pocket = str(i + 1)

            pocket_id = df.at[index, "cav_id"]

            pocket_chainid_lst = lst_unique(
                (df.at[index, "name_chain_1"], df.at[index, "name_chain_2"])
            )

            pocket_chainid = lst_to_str(
                pocket_chainid_lst,
                join_txt="_",
            )

            cont_path = get_file_path(
                f"pocket{pocket_id}_atm.pdb", dir_path=run_dir_path
            )
            sphere_path = get_file_path(
                f"pocket{pocket_id}_vert.pqr", dir_path=run_dir_path
            )

            sphere_obj = f"sphere_{pocket}"
            cmd.load(sphere_path, sphere_obj)

            cmd.alter(sphere_obj, f"chain='{pocket_chainid}'")

            pocket_obj = f"pocket_{pocket}"
            cmd.create(pocket_obj, f"{obj} {sphere_obj}")

            pocket_path = get_pocket_path(
                obj,
                pocket,
                dir_path=pocket_dir,
            )
            cmd.save(pocket_path, pocket_obj)

            structure = load_coord(cont_path)

            pocket_residue_lst = get_residues(structure)

            if chainid_lst is not None:
                pocket_residue_lst = [
                    x for x in pocket_residue_lst if get_reschainid(x) in chainid_lst
                ]

            pocket_cont_lst = sort_lst(
                lst_unique(
                    [
                        get_resnum(x)
                        for x in pocket_residue_lst
                        if (get_resnum(x) < 50000) and is_aa(x)
                    ]
                )
            )

            pocket_lig = df.at[index, "lig_het_tag"]

            pocket_status = pocket_unbound_name
            pocket_type = pocket_unbound_name
            if pd.isna(pocket_lig):
                pocket_lig = "STP"
            else:
                pocket_status = pocket_bound_name
                is_pharm = True
                for lig_col in lig_col_lst:
                    if pocket_lig in lig_lst_dict[lig_col]:
                        pocket_type = lig_col
                        is_pharm = False
                if is_pharm:
                    pocket_type = pharm_lig_col
                pocket_type = title_str(pocket_type)

            pocket_dict[coord_path][pocket] = dict()

            pocket_dict[coord_path][pocket][pocket_volume_col] = df.at[index, "volume"]
            pocket_dict[coord_path][pocket][pocket_score_col] = df.at[
                index, "drug_score"
            ]
            pocket_dict[coord_path][pocket][pocket_status_col] = pocket_status
            pocket_dict[coord_path][pocket][pocket_type_col] = pocket_type
            pocket_dict[coord_path][pocket][pocket_lig_col] = pocket_lig
            pocket_dict[coord_path][pocket][pocket_cont_col] = lst_to_str(
                pocket_cont_lst
            )
            pocket_dict[coord_path][pocket][pocket_path_col] = pocket_path

    delete_path(run_path)
    delete_path(run_dir_path)

    return pocket_dict


def prep_pocket(
    coord_paths,
    pocket_dir=None,
    pocket_json_path=None,
    chainid_dict=None,
    update_pocket=False,
    num_cpu=1,
):

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

    if update_pocket:
        delete_path(pocket_dir_path)
    append_path(pocket_dir_path)

    coord_path_lst = type_lst(coord_paths)

    pocket_dict = dict()

    if num_cpu == 1:
        for coord_path in tqdm(
            coord_path_lst,
            desc="Preparing pockets",
            position=0,
            leave=True,
        ):
            pocket_dict = merge_dicts(
                [
                    pocket_dict,
                    run_fpocket(
                        coord_path,
                        pocket_dir=pocket_dir,
                        chainid_lst=chainid_dict[coord_path],
                    ),
                ]
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    run_fpocket,
                    coord_path,
                    pocket_dir=pocket_dir,
                    chainid_lst=chainid_dict[coord_path],
                )
                for coord_path in coord_path_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Preparing pockets",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                pocket_dict = merge_dicts([pocket_dict, job.result()])

    print("Prepared pockets!")

    if pocket_json_path is not None:
        save_json(pocket_json_path, pocket_dict)
    else:
        return pocket_dict