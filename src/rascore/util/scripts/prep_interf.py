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
import pymol2

from ..functions.table import merge_dicts
from ..functions.path import (
    get_file_name,
    get_interf_path,
    get_dir_path,
    delete_path,
    append_path,
    save_json,
    interf_str,
)
from ..functions.lst import (
    lst_to_str,
    sort_lst,
    lst_unique,
    type_lst,
    calc_simpson,
    calc_jaccard,
)
from ..functions.col import (
    interf_path_col,
    bound_interf_chainid_col,
    interf_area_col,
    cb_cont_col,
    atomid_cont_col,
    total_cb_cont_col,
    total_atomid_cont_col,
    iso_col,
)


byres_str = "byres"
within_str = "within"
of_str = "of"
cb_str = "and name CB"

max_cb_dist = 12.0
max_atomid_dist = 5.0


def save_interf(
    coord_path,
    interf_dir=None,
    sym_dist=5.0,
    min_area=200,
    iso_min_simi=0.7,
    use_simpson=True,
    chainid_lst=None,
):

    file_name = get_file_name(coord_path)

    if ".cif" in file_name:
        file_format = "cif"
        return_pdb = False
    elif ".pdb" in file_name:
        file_format = "pdb"
        return_pdb = True

    obj = file_name.replace(f".{file_format}", "")

    if ".gz" in file_name:
        obj = obj.replace(".gz", "")

    with pymol2.PyMOL() as pymol:

        cmd = pymol.cmd

        cmd.load(coord_path, obj)

        cmd.set("dot_solvent", 1)

        if chainid_lst is None:
            chainid_lst = cmd.get_chains(obj)
        else:
            cmd.remove(f"{obj} and not chain {lst_to_str(chainid_lst,join_txt='+')}")

        cmd.split_chains(obj)

        for chainid in chainid_lst:
            chainid_obj = f"{obj}_{chainid}"

            cmd.symexp(f"{chainid_obj}_", chainid_obj, chainid_obj, sym_dist)

        sym_obj_lst = cmd.get_object_list("all")
        sym_obj_lst.remove(obj)

        interf_dict = {coord_path: dict()}

        for chainid in chainid_lst:
            chainid_obj = f"{obj}_{chainid}"

            chainid_area = cmd.get_area(f"{chainid_obj} and polymer")

            interf = 1

            interf_dict[coord_path][chainid] = dict()

            past_cb_lst = list()
            past_atomid_lst = list()

            for sym_obj in sym_obj_lst:
                if chainid_obj != sym_obj:
                    return_chainid = cmd.get_chains(sym_obj)[0]

                    sym_chainid = f"{return_chainid}_{interf}"

                    cmd.alter(sym_obj, f"chain='{sym_chainid}'")

                    interf_obj = f"{chainid_obj}_{interf}"

                    cmd.create(interf_obj, f"{chainid_obj} {sym_obj}")

                    cmd.alter(sym_obj, f"chain='{return_chainid}'")

                    sym_area = cmd.get_area(f"{sym_obj} and polymer")
                    complex_area = cmd.get_area(f"{interf_obj} and polymer")

                    interf_area = ((chainid_area + sym_area) - complex_area) / 2

                    if interf_area >= min_area:

                        chainid_1_sele = f"{interf_obj} and chain {chainid} and polymer"
                        chainid_2_sele = (
                            f"{interf_obj} and chain {sym_chainid} and polymer"
                        )

                        chainid_cb_iterate = lst_to_str(
                            [
                                byres_str,
                                chainid_1_sele,
                                cb_str,
                                within_str,
                                max_cb_dist,
                                of_str,
                                chainid_2_sele,
                                cb_str,
                            ],
                            join_txt=" ",
                        )
                        chainid_atomid_iterate = lst_to_str(
                            [
                                byres_str,
                                chainid_1_sele,
                                within_str,
                                max_atomid_dist,
                                of_str,
                                chainid_2_sele,
                            ],
                            join_txt=" ",
                        )

                        sym_cb_iterate = lst_to_str(
                            [
                                byres_str,
                                chainid_2_sele,
                                cb_str,
                                within_str,
                                max_cb_dist,
                                of_str,
                                chainid_1_sele,
                                cb_str,
                            ],
                            join_txt=" ",
                        )

                        sym_atomid_iterate = lst_to_str(
                            [
                                byres_str,
                                chainid_2_sele,
                                within_str,
                                max_atomid_dist,
                                of_str,
                                chainid_1_sele,
                            ],
                            join_txt=" ",
                        )

                        cont_dict = {
                            "chainid_cb_lst": list(),
                            "chainid_atomid_lst": list(),
                            "sym_cb_lst": list(),
                            "sym_atomid_lst": list(),
                        }

                        cmd.iterate(
                            chainid_cb_iterate,
                            "chainid_cb_lst.append(resi)",
                            space=cont_dict,
                        )
                        cmd.iterate(
                            chainid_atomid_iterate,
                            "chainid_atomid_lst.append(resi)",
                            space=cont_dict,
                        )

                        cmd.iterate(
                            sym_cb_iterate,
                            "sym_cb_lst.append(resi)",
                            space=cont_dict,
                        )
                        cmd.iterate(
                            sym_atomid_iterate,
                            "sym_atomid_lst.append(resi)",
                            space=cont_dict,
                        )

                        chainid_cb_lst = sort_lst(
                            lst_unique(cont_dict["chainid_cb_lst"])
                        )

                        chainid_atomid_lst = sort_lst(
                            lst_unique(cont_dict["chainid_atomid_lst"])
                        )

                        sym_cb_lst = sort_lst(lst_unique(cont_dict["sym_cb_lst"]))
                        sym_atomid_lst = sort_lst(
                            lst_unique(cont_dict["sym_atomid_lst"])
                        )

                        cb_lst = sort_lst(lst_unique(chainid_cb_lst + sym_cb_lst))
                        atomid_lst = sort_lst(
                            lst_unique(chainid_atomid_lst + sym_atomid_lst)
                        )

                        if not (cb_lst in past_cb_lst or atomid_lst in past_atomid_lst):

                            past_cb_lst.append(cb_lst)
                            past_atomid_lst.append(atomid_lst)

                            total_cb_cont = len(chainid_cb_lst)
                            total_atomid_cont = len(chainid_atomid_lst)

                            if (
                                total_cb_cont >= 10 and total_atomid_cont >= 1
                            ) or total_atomid_cont >= 5:

                                interf_path = get_interf_path(
                                    obj,
                                    chainid,
                                    interf,
                                    dir_path=interf_dir,
                                    return_pdb=return_pdb,
                                )
                                cmd.save(interf_path, interf_obj)

                                iso_status = False
                                if use_simpson:
                                    cb_simi = calc_simpson(chainid_cb_lst, sym_cb_lst)
                                    atomid_simi = calc_simpson(
                                        chainid_atomid_lst, sym_atomid_lst
                                    )
                                else:
                                    cb_simi = calc_jaccard(chainid_cb_lst, sym_cb_lst)
                                    atomid_simi = calc_jaccard(
                                        chainid_atomid_lst, sym_atomid_lst
                                    )
                                if (
                                    cb_simi >= iso_min_simi
                                    and atomid_simi >= iso_min_simi
                                ):
                                    iso_status = True

                                interf_dict[coord_path][chainid][interf] = dict()

                                interf_dict[coord_path][chainid][interf][
                                    interf_path_col
                                ] = interf_path
                                interf_dict[coord_path][chainid][interf][
                                    bound_interf_chainid_col
                                ] = sym_chainid
                                interf_dict[coord_path][chainid][interf][
                                    interf_area_col
                                ] = interf_area
                                interf_dict[coord_path][chainid][interf][
                                    cb_cont_col
                                ] = lst_to_str(chainid_cb_lst)
                                interf_dict[coord_path][chainid][interf][
                                    atomid_cont_col
                                ] = lst_to_str(chainid_atomid_lst)
                                interf_dict[coord_path][chainid][interf][
                                    total_cb_cont_col
                                ] = total_cb_cont
                                interf_dict[coord_path][chainid][interf][
                                    total_atomid_cont_col
                                ] = total_atomid_cont
                                interf_dict[coord_path][chainid][interf][
                                    iso_col
                                ] = iso_status

                                interf += 1

    return interf_dict


def prep_interf(
    coord_paths,
    interf_dir=None,
    interf_json_path=None,
    sym_dist=5.0,
    min_area=200,
    iso_min_simi=0.7,
    use_simpson=True,
    chainid_dict=None,
    update_interf=False,
    num_cpu=1,
):

    interf_path = get_dir_path(dir_str=interf_str, dir_path=interf_dir)

    if update_interf:
        delete_path(interf_path)
    append_path(interf_path)

    coord_path_lst = type_lst(coord_paths)

    interf_dict = dict()

    if num_cpu == 1:
        for coord_path in tqdm(
            coord_path_lst,
            desc="Preparing interfaces",
            position=0,
            leave=True,
        ):
            interf_dict = merge_dicts(
                [
                    interf_dict,
                    save_interf(
                        coord_path,
                        interf_dir=interf_dir,
                        sym_dist=sym_dist,
                        min_area=min_area,
                        chainid_lst=chainid_dict[coord_path],
                    ),
                ]
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    save_interf,
                    coord_path,
                    interf_dir=interf_dir,
                    sym_dist=sym_dist,
                    min_area=min_area,
                    iso_min_simi=iso_min_simi,
                    use_simpson=use_simpson,
                    chainid_lst=chainid_dict[coord_path],
                )
                for coord_path in coord_path_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Preparing interfaces",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                interf_dict = merge_dicts([interf_dict, job.result()])

    print("Prepared interfaces!")

    if interf_json_path is not None:
        save_json(interf_json_path, interf_dict)
    else:
        return interf_dict
