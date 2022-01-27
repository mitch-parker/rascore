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
    calc_simpson,
    calc_jaccard,
)
from ..functions.table import mask_equal, lst_col, title_str, merge_dicts
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
from ..functions.file import pocket_table_file
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


def run_dpocket(pharm_dict, chainid_dict=None, pocket_dir=None):

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

    dp_file_path = get_file_path(
        pocket_table_file, dir_str=pocket_str, dir_path=pocket_dir
    )

    coord_path_dict = dict()
    run_chainid_dict = dict()

    with open(dp_file_path, "w") as file:
        for coord_path in list(pharm_dict.keys()):

            coord_name = get_file_name(coord_path)
            coord_name = coord_name.replace("_core", "")

            run_path = get_file_path(
                coord_name, dir_str=pocket_str, dir_path=pocket_dir
            )

            coord_path_dict[run_path] = coord_path
            if chainid_dict is not None:
                run_chainid_dict[run_path] = chainid_dict[coord_path]

            copy_path(coord_path, run_path)

            file.write(f"{run_path}\t{pharm_dict[coord_path]}\n")

    os.system(f"dpocket -f {dp_file_path} -o {pocket_dir_path}/pocket")

    fp_file_path = get_file_path("fp.txt", dir_str=pocket_str, dir_path=pocket_dir)
    fpn_file_path = get_file_path("fpn.txt", dir_str=pocket_str, dir_path=pocket_dir)
    exp_file_path = get_file_path("exp.txt", dir_str=pocket_str, dir_path=pocket_dir)

    df = pd.read_csv(
        fp_file_path,
        delim_whitespace=True,
        dtype=str,
    )

    pocket_dict = dict()

    for run_path in tqdm(
        lst_col(df, "pdb", unique=True),
        desc="Preparing pockets",
        position=0,
        leave=True,
    ):

        run_df = mask_equal(df, "pdb", run_path)

        if len(run_df) > 1:
            run_df = mask_equal(run_df, "overlap", run_df["overlap"].max())

        run_name = get_file_name(run_path)
        obj = run_name.replace(".pdb", "")

        run_dir_path = f"{pocket_dir_path}/{obj}_out"
        run_pockets_dir_path = f"{run_dir_path}/pockets"

        info_file_path = get_file_path(f"{obj}_info.txt", dir_path=run_dir_path)

        add_pocket = False
        with open(info_file_path, "r") as file:
            for line in file.readlines():
                if line[:6] == "Pocket":
                    pocket_id = line[7:8]
                    count = 0
                for info, col in pocket_info_dict.items():
                    info_str = f"{info} :"
                    if info_str in line:
                        info_val = round(
                            float(line.split(info_str)[1].replace(" ", ""))
                        )
                        col_val = round(float(run_df.at[0, col]))
                        if col_val == info_val:
                            count += 1
                if count == 3:
                    add_pocket = True
                    break

        if add_pocket:
            cont_path = get_file_path(
                f"pocket{pocket_id}_atm.pdb", dir_path=run_pockets_dir_path
            )
            structure = load_coord(cont_path)

            pocket_residue_lst = get_residues(structure)

            if chainid_dict is not None:
                chainid_lst = run_chainid_dict[run_path]
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

            pocket_lig = run_df.at[0, "lig"]
            pocket_path = get_pocket_path(
                obj.replace("pocket_", ""),
                pocket_lig,
                dir_path=pocket_dir,
            )

            coord_path = coord_path_dict[run_path]

            copy_path(coord_path, pocket_path)

            pocket_dict[coord_path] = {pocket_lig: dict()}

            pocket_dict[coord_path][pocket_lig][pocket_volume_col] = run_df.at[
                0, "pock_vol"
            ]
            pocket_dict[coord_path][pocket_lig][pocket_score_col] = run_df.at[
                0, "drug_score"
            ]
            pocket_dict[coord_path][pocket_lig][pocket_status_col] = pocket_bound_name
            pocket_dict[coord_path][pocket_lig][pocket_type_col] = title_str(
                pharm_lig_col
            )
            pocket_dict[coord_path][pocket_lig][pocket_lig_col] = pocket_lig
            pocket_dict[coord_path][pocket_lig][pocket_cont_col] = lst_to_str(
                pocket_cont_lst
            )
            pocket_dict[coord_path][pocket_lig][pocket_path_col] = pocket_path

        delete_path(run_path)
        delete_path(run_dir_path)

    delete_path(dp_file_path)
    delete_path(fp_file_path)
    delete_path(fpn_file_path)
    delete_path(exp_file_path)

    return pocket_dict


def run_fpocket(
    coord_path,
    pocket_dir=None,
    chainid_lst=None,
    check_lig=False,
    check_dist=3.5,
    check_min_simi=0.1,
    use_simpson=True,
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

            pocket = i + 1

            pocket_id = df.at[index, "cav_id"]

            pocket_chainid_lst = lst_unique(
                (df.at[index, "name_chain_1"], df.at[index, "name_chain_2"])
            )

            if check_lig:
                pocket_chainid_sele = lst_to_str(pocket_chainid_lst, join_txt="+")

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
                type_pocket = True

                if check_lig:
                    lig_cont_str = f"byres polymer and chain {pocket_chainid_sele} within {check_dist} of resname {pocket_lig}"
                    lig_cont_dict = {
                        "lig_cont_lst": list(),
                    }

                    cmd.iterate(
                        lig_cont_str,
                        "lig_cont_lst.append(resi)",
                        space=lig_cont_dict,
                    )

                    lig_cont_lst = [
                        int(x) for x in lst_unique(lig_cont_dict["lig_cont_lst"])
                    ]

                    if use_simpson:
                        check_simi = calc_simpson(pocket_cont_lst, lig_cont_lst)
                    else:
                        check_simi = calc_jaccard(pocket_cont_lst, lig_cont_lst)

                    if check_simi <= check_min_simi:
                        pocket_lig = "STP"
                        type_pocket = False

                if type_pocket:
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
    pharm_dict=None,
    chainid_dict=None,
    check_lig=False,
    check_dist=3.5,
    check_min_simi=0.1,
    use_simpson=True,
    update_pocket=False,
    num_cpu=1,
):

    pocket_dir_path = get_dir_path(dir_str=pocket_str, dir_path=pocket_dir)

    if update_pocket:
        delete_path(pocket_dir_path)
    append_path(pocket_dir_path)

    coord_path_lst = type_lst(coord_paths)

    pocket_dict = dict()

    if pharm_dict is not None:
        coord_path_lst = [x for x in coord_path_lst if x not in list(pharm_dict.keys())]
        pocket_dict = merge_dicts(
            [
                pocket_dict,
                run_dpocket(
                    pharm_dict,
                    chainid_dict=chainid_dict,
                    pocket_dir=pocket_dir,
                ),
            ]
        )

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
                        check_lig=check_lig,
                        check_dist=check_dist,
                        check_min_simi=check_min_simi,
                        use_simpson=use_simpson,
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
                    check_lig=check_lig,
                    check_dist=check_dist,
                    check_min_simi=check_min_simi,
                    use_simpson=use_simpson,
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