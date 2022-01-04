# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import subprocess
import io

import pandas as pd

from tqdm import tqdm
import concurrent.futures
import pymol2

from functions import *


def run_fpocket(
    coord_path,
    pocket_dir=None,
    chainid_lst=None,
):

    coord_name = get_file_name(coord_path)
    coord_name = coord_name.replace("_core", "")

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

    run_path = get_file_path(coord_name, dir_str=pocket_str, dir_path=pocket_dir)

    copy_path(coord_path, run_path)

    cmd_lst = ["fpocket", "-f", run_path, "-d"]

    df = pd.read_csv(
        io.StringIO(
            subprocess.Popen(cmd_lst, stdout=subprocess.PIPE)
            .communicate()[0]
            .decode("utf-8")
        ),
        sep=" ",
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

            pocket = i + 1

            pocket_id = df.at[index, "cav_id"]

            pocket_chainid = lst_to_str(
                lst_unique(
                    (df.at[index, "name_chain_1"], df.at[index, "name_chain_2"])
                ),
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

            pocket_dict[coord_path][pocket] = dict()

            pocket_dict[coord_path][pocket][pocket_volume_col] = df.at[index, "volume"]
            pocket_dict[coord_path][pocket][pocket_score_col] = df.at[
                index, "drug_score"
            ]
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
    num_cpu=1,
):

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

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