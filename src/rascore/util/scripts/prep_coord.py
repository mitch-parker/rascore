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
import pandas as pd
import xml.etree.ElementTree as ET
import gzip
from tqdm import tqdm
import concurrent.futures
from Bio.PDB import Select
import pymol2

from ..functions.coord import (
    load_coord,
    load_cif_dict,
    save_coord,
    build_pdb_dict,
    build_pdb_code_lst,
    get_reschainid,
    get_resid,
    get_atomid,
    get_modelid,
    get_neighbors,
    get_chainid,
    get_chain_cont,
    get_resnum,
    get_resname,
    is_aa,
    is_het,
    is_wat,
    join_seq_lst,
    get_seq_lst,
    get_residue_cont,
)
from ..functions.path import (
    path_exists,
    get_renum_path,
    get_dir_path,
    save_table,
    save_json,
    get_core_path,
    get_rcsb_path,
    get_sifts_path,
    search_dir,
    append_path,
    delete_path,
    core_str,
    sifts_str,
    rcsb_str,
    renum_str,
    rcsb_assembly_str,
    renum_assembly_str,
)
from ..functions.lst import lst_to_str, res_to_str, type_lst, lst_unique, res_to_lst
from ..functions.table import get_str_num, merge_dicts, merge_tables
from ..functions.col import (
    modelid_col,
    bound_lig_col,
    bound_prot_chainid_col,
    range_col,
    pdb_code_col,
    chainid_col,
    rcsb_path_col,
    renum_path_col,
    sifts_path_col,
    rcsb_assembly_path_col,
    renum_assembly_path_col,
    core_path_col,
    pdb_id_col,
)

max_cb_dist = 12.0
max_atomid_dist = 5.0


class ChainSelect(Select):
    def __init__(self, sele_dict):
        self.sele_dict = sele_dict

    def accept_residue(self, residue):

        return self.sele_dict[get_reschainid(residue)][get_resid(residue)]

    def accept_atom(self, atom):

        if get_atomid(atom)[0] == "H":
            return 0
        else:
            return 1


def run_pdb_renum(
    pdb_code_lst,
    renum_script_path,
    rcsb_dir=None,
    sifts_dir=None,
    renum_dir=None,
    num_cpu=1,
):

    renum_pdb_code_lst = [
        x
        for x in pdb_code_lst
        if not path_exists(get_renum_path(x, dir_path=renum_dir))
    ]

    if len(renum_pdb_code_lst) > 0:
        pdb_code_str = lst_to_str(
            renum_pdb_code_lst,
            join_txt=" ",
        )

        cmd_lst = [
            f"python {renum_script_path}",
            f"-rfla {pdb_code_str}",
            "-mmCIF",
            f"-sipm {get_dir_path(dir_str=rcsb_str, dir_path=rcsb_dir)}",
            f"-sips {get_dir_path(dir_str=sifts_str,dir_path=sifts_dir)}",
            f"-sopm {get_dir_path(dir_str=renum_str, dir_path=renum_dir)}",
            "-offz",
            f"-nproc {num_cpu}",
        ]

        cmd_str = lst_to_str(cmd_lst, join_txt=" ")

        os.system(cmd_str)

        assembly_cmd_lst = [
            f"python {renum_script_path}",
            f"-rfla {pdb_code_str}",
            "-mmCIF_assembly",
            f"-sipma {get_dir_path(dir_str=rcsb_assembly_str, dir_path=rcsb_dir)}",
            f"-sips {get_dir_path(dir_str=sifts_str,dir_path=sifts_dir)}",
            f"-sopma {get_dir_path(dir_str=renum_assembly_str, dir_path=renum_dir)}",
            "-offz",
            f"-nproc {num_cpu}",
        ]

        assembly_cmd_str = lst_to_str(assembly_cmd_lst, join_txt=" ")

        os.system(assembly_cmd_str)


def get_sifts_dict(sifts_path):

    tree = ET.parse(gzip.open(sifts_path, "rt"))
    root = tree.getroot()
    base = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"

    sifts_dict = dict()
    uniprot_dict = dict()

    pdb_chainid_lst = list()
    uniprot_chainid_lst = list()

    for entity in root:
        if entity.tag == (str(base) + "entity"):
            for segment in entity:
                if segment.tag == (str(base) + "segment"):
                    for listResidue in segment:
                        if listResidue.tag == (str(base) + "listResidue"):
                            for residue in listResidue:
                                if residue.tag == (str(base) + "residue"):
                                    if residue.get("dbSource") == "PDBe":

                                        pdbe_resid = residue.get("dbResNum")

                                        pdb_resid = "null"
                                        uniprot_resid = "null"

                                        no_uniprot = True

                                        for crossRefDb in residue:
                                            if crossRefDb.get("dbSource") == "PDB":
                                                pdb_code = crossRefDb.get(
                                                    "dbAccessionId"
                                                )
                                                chainid = crossRefDb.get("dbChainId")

                                                if chainid not in pdb_chainid_lst:
                                                    pdb_chainid_lst.append(chainid)

                                                pdb_resid = crossRefDb.get("dbResNum")

                                            if crossRefDb.get("dbSource") == "UniProt":

                                                if chainid not in uniprot_chainid_lst:
                                                    uniprot_chainid_lst.append(chainid)

                                                no_uniprot = False
                                                uniprot_resid = crossRefDb.get(
                                                    "dbResNum"
                                                )

                                        if no_uniprot:
                                            uniprot_resid = str(
                                                get_str_num(pdbe_resid) + 50000
                                            )

                                            if chainid not in list(uniprot_dict.keys()):
                                                uniprot_dict[chainid] = dict()

                                            uniprot_dict[chainid][
                                                uniprot_resid
                                            ] = pdb_resid

                                        if (
                                            pdb_resid != "null"
                                            and uniprot_resid != "null"
                                        ):
                                            if pdb_code not in list(sifts_dict.keys()):
                                                sifts_dict[pdb_code] = dict()

                                            if chainid not in list(
                                                sifts_dict[pdb_code].keys()
                                            ):
                                                sifts_dict[pdb_code][chainid] = dict()

                                            sifts_dict[pdb_code][chainid][
                                                pdb_resid
                                            ] = uniprot_resid

    for chainid in pdb_chainid_lst:
        if chainid not in uniprot_chainid_lst:
            if chainid in list(uniprot_dict.keys()):
                for pdb_resid in list(sifts_dict[pdb_code][chainid].keys()):
                    uniprot_resid = sifts_dict[pdb_code][chainid][pdb_resid]
                    if int(uniprot_resid) > 50000:
                        sifts_dict[pdb_code][chainid][pdb_resid] = uniprot_dict[
                            chainid
                        ][uniprot_resid]

    return sifts_dict


def build_sifts_map(sifts_path_lst, sifts_json_path, num_cpu=1):

    sifts_dict = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:

        job_lst = [
            executor.submit(get_sifts_dict, sifts_path) for sifts_path in sifts_path_lst
        ]

        for job in tqdm(
            concurrent.futures.as_completed(job_lst),
            desc="Preparing SIFTS numbering",
            total=len(job_lst),
            miniters=1,
            position=0,
            leave=True,
        ):

            sifts_dict = merge_dicts([sifts_dict, job.result()])

    save_json(sifts_json_path, sifts_dict)


def build_sele_dict(
    structure,
    curr_chainid,
    chainid_lst,
    max_lig_dist=4,
    min_lig_cont=5,
    prot_resid_lst=None,
    lig_resid_lst=None,
    bound_lig_lst=None,
    bound_chainid_lst=None,
):

    sele_dict = dict()
    sele_df = pd.DataFrame()

    for model in structure:

        modelid = get_modelid(model)

        sele_dict[modelid] = dict()

        curr_chain = model[curr_chainid]

        neighbors = get_neighbors(curr_chain)

        resnum_lst = list()
        bound_residue_dict = dict()
        bound_chain_dict = dict()

        for chain in model:

            chainid = get_chainid(chain)

            sele_dict[modelid][chainid] = dict()

            if chainid not in chainid_lst:

                sele = 0

                check_chainid = True
                if bound_chainid_lst is not None:
                    if chainid not in bound_chainid_lst:
                        check_chainid = False

                if check_chainid:

                    atomid_chain_cont = [
                        x
                        for x in get_chain_cont(
                            neighbors,
                            chain,
                            max_dist=max_atomid_dist,
                        )
                        if is_aa(x.get_parent())
                        and get_reschainid(x.get_parent()) == curr_chainid
                        and get_resnum(x.get_parent()) < 50000
                    ]

                    cb_chain_cont = [
                        x
                        for x in get_chain_cont(
                            neighbors,
                            chain,
                            max_dist=max_cb_dist,
                        )
                        if is_aa(x.get_parent())
                        and get_reschainid(x.get_parent()) == curr_chainid
                        and get_resnum(x.get_parent()) < 50000
                        and get_atomid(x) == "CB"
                    ]

                    if prot_resid_lst is not None:
                        atomid_chain_cont = [
                            x
                            for x in atomid_chain_cont
                            if get_resnum(x.get_parent()) in prot_resid_lst
                        ]
                        cb_chain_cont = [
                            x
                            for x in cb_chain_cont
                            if get_resnum(x.get_parent()) in prot_resid_lst
                        ]

                    atomid_chain_cont = len(atomid_chain_cont)
                    cb_chain_cont = len(cb_chain_cont)

                    chain_seq = join_seq_lst(get_seq_lst(chain))

                    if (
                        cb_chain_cont >= 10 and atomid_chain_cont >= 1
                    ) or atomid_chain_cont >= 5:

                        if chain_seq not in list(bound_chain_dict.keys()):
                            bound_chain_dict[chain_seq] = dict()
                        bound_chain_dict[chain_seq][chainid] = cb_chain_cont

                for residue in chain:
                    resid = get_resid(residue)
                    sele_dict[modelid][chainid][resid] = sele
            else:
                for residue in chain:

                    resid = get_resid(residue)
                    sele = 0

                    if is_aa(residue):

                        if chainid == curr_chainid:
                            resnum = get_resnum(residue)
                            if resnum < 50000:
                                resnum_lst.append(resnum)
                            sele = 1

                    elif is_het(residue):

                        atomid_residue_cont = [
                            x
                            for x in get_residue_cont(
                                neighbors,
                                residue,
                                max_dist=max_lig_dist,
                            )
                            if (
                                is_het(x.get_parent())
                                and get_reschainid(x.get_parent()) == curr_chainid
                            )
                            or (
                                is_aa(x.get_parent())
                                and get_reschainid(x.get_parent()) == curr_chainid
                                and get_resnum(x.get_parent()) < 50000
                            )
                        ]

                        if lig_resid_lst is not None:
                            atomid_residue_cont = [
                                x
                                for x in atomid_residue_cont
                                if is_het(x.get_parent())
                                or (
                                    is_aa(x.get_parent())
                                    and get_resnum(x.get_parent()) in lig_resid_lst
                                )
                            ]

                        atomid_residue_cont = len(atomid_residue_cont)

                        if atomid_residue_cont >= min_lig_cont:
                            resname = get_resname(residue)
                            if resname not in list(bound_residue_dict.keys()):
                                bound_residue_dict[resname] = dict()
                            bound_residue_dict[resname][
                                (chainid, resid)
                            ] = atomid_residue_cont

                    elif is_wat(residue):
                        if chainid == curr_chainid:
                            sele = 1

                    sele_dict[modelid][chainid][resid] = sele

        bound_prot_chainid_lst = list()
        for chain_seq in list(bound_chain_dict.keys()):
            bound_chain_lst = list(bound_chain_dict[chain_seq].keys())
            bound_chain_cont_lst = [
                bound_chain_dict[chain_seq][x] for x in bound_chain_lst
            ]
            bound_chain = bound_chain_lst[
                bound_chain_cont_lst.index(max(bound_chain_cont_lst))
            ]
            bound_prot_chainid_lst.append(bound_chain)
            for residue in model[bound_chain]:
                resid = get_resid(residue)
                sele_dict[modelid][bound_chain][resid] = 1

        bound_lig_lst = list()
        for resname in list(bound_residue_dict.keys()):
            bound_residue_lst = list(bound_residue_dict[resname].keys())
            bound_residue_cont_lst = [
                bound_residue_dict[resname][x] for x in bound_residue_lst
            ]
            bound_residue = bound_residue_lst[
                bound_residue_cont_lst.index(max(bound_residue_cont_lst))
            ]
            bound_lig_lst.append(resname)
            sele_dict[modelid][bound_residue[0]][bound_residue[1]] = 1

        sele_df.at[modelid, modelid_col] = modelid
        sele_df.at[modelid, bound_lig_col] = lst_to_str(bound_lig_lst)
        sele_df.at[modelid, bound_prot_chainid_col] = lst_to_str(bound_prot_chainid_lst)
        sele_df.at[modelid, range_col] = res_to_str(resnum_lst)

    return sele_dict, sele_df


def isolate_chains(
    pdb_code,
    pdb_dict,
    core_dir=None,
    rcsb_dir=None,
    sifts_dir=None,
    renum_dir=None,
    max_lig_dist=4,
    min_lig_cont=5,
    prot_resid_lst=None,
    lig_resid_lst=None,
    add_h=False,
    add_h_his=False,
    bound_chainid_dict=None,
    all_models=False,
    update_coords=False,
):

    rcsb_path = get_rcsb_path(pdb_code, dir_path=rcsb_dir)
    renum_path = get_renum_path(pdb_code, dir_path=renum_dir)
    sifts_path = get_sifts_path(pdb_code, dir_path=sifts_dir)

    chainid_lst = pdb_dict[pdb_code]

    structure = load_coord(renum_path)

    if not all_models:
        structure = type_lst(structure[0])

    chain_df = pd.DataFrame()

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd

        for chainid in chainid_lst:

            bound_chainid_lst = None
            if bound_chainid_dict is not None:
                pdb_id = f"{pdb_code}{chainid}"
                if pdb_id in (bound_chainid_dict.keys()):
                    bound_chainid_lst = bound_chainid_dict[pdb_id]

            rcsb_assembly_path = "None"
            renum_assembly_path = "None"

            if bound_chainid_lst is None:

                renum_assembly_dir = get_dir_path(
                    dir_str=renum_assembly_str, dir_path=renum_dir
                )
                rcsb_assembly_dir = get_dir_path(
                    dir_str=rcsb_assembly_str, dir_path=rcsb_dir
                )

                for renum_assembly_file in search_dir(renum_assembly_dir, pdb_code):

                    renum_assembly_path = f"{renum_assembly_dir}/{renum_assembly_file}"
                    rcsb_assembly_path = f"{rcsb_assembly_dir}/{renum_assembly_file.replace('_renum.cif','.cif.gz')}"

                    assembly_dict = load_cif_dict(renum_assembly_path)

                    assembly_chainid_lst = lst_unique(
                        assembly_dict["_pdbe_chain_remapping.orig_auth_asym_id"]
                    )

                    if chainid in assembly_chainid_lst:
                        if (
                            len(
                                [
                                    x
                                    for x in assembly_chainid_lst
                                    if x in chainid_lst and x != chainid
                                ]
                            )
                            == 0
                        ):
                            bound_chainid_lst = [
                                x for x in assembly_chainid_lst if x != chainid
                            ]

            sele_dict, sele_df = build_sele_dict(
                structure,
                chainid,
                chainid_lst,
                max_lig_dist=max_lig_dist,
                min_lig_cont=min_lig_cont,
                prot_resid_lst=prot_resid_lst,
                lig_resid_lst=lig_resid_lst,
                bound_chainid_lst=bound_chainid_lst,
            )

            sele_df[pdb_code_col] = pdb_code
            sele_df[chainid_col] = chainid
            sele_df[rcsb_path_col] = rcsb_path
            sele_df[renum_path_col] = renum_path
            sele_df[sifts_path_col] = sifts_path
            sele_df[rcsb_assembly_path_col] = rcsb_assembly_path
            sele_df[renum_assembly_path_col] = renum_assembly_path

            for modelid in list(sele_dict.keys()):

                modelid_str = modelid
                if not all_models:
                    modelid_str = None

                cif_path = get_core_path(
                    pdb_code, chainid, modelid=modelid_str, dir_path=core_dir
                )
                pdb_path = get_core_path(
                    pdb_code,
                    chainid,
                    modelid=modelid_str,
                    return_pdb=True,
                    dir_path=core_dir,
                )

                update_cif = True
                if not update_coords:
                    if path_exists(cif_path):
                        update_cif = False

                update_pdb = True
                if not update_coords:
                    if path_exists(cif_path):
                        update_pdb = False

                sele_df.at[int(modelid), core_path_col] = cif_path

                if update_cif:
                    chain_sele = ChainSelect(sele_dict[modelid])
                    save_coord(cif_path, structure[int(modelid)], sele=chain_sele)

                if update_pdb:
                    pymol_obj = f"{pdb_code}{chainid}"

                    if all_models:
                        pymol_obj += str(modelid)

                    cmd.load(cif_path, pymol_obj)
                    cmd.save(pdb_path, pymol_obj)

                if add_h:
                    pdb_h_path = get_core_path(
                        pdb_code,
                        chainid,
                        modelid=modelid_str,
                        return_pdb=True,
                        add_h=True,
                        dir_path=core_dir,
                    )

                    update_pdb_h = True
                    if not update_coords:
                        if path_exists(pdb_h_path):
                            update_pdb_h = False

                    if update_pdb_h:
                        cmd_lst = ["reduce", pdb_path, ">", pdb_h_path, "-NOHETh"]

                        if add_h_his:
                            cmd_lst.append("-HIS")

                        os.system(lst_to_str(cmd_lst, join_txt=" "))

            sele_df[pdb_id_col] = sele_df[pdb_code_col].map(str) + sele_df[
                chainid_col
            ].map(str)

            if all_models:
                sele_df[pdb_id_col] += sele_df[modelid_col].map(str)

            sele_df[modelid_col] = 0

            chain_df = pd.concat([chain_df, sele_df], sort=False)

    return chain_df


def run_pdb_chain(
    pdb_ids,
    core_dir=None,
    rcsb_dir=None,
    sifts_dir=None,
    renum_dir=None,
    max_lig_dist=4,
    min_lig_cont=5,
    prot_resids=None,
    lig_resids=None,
    add_h=False,
    add_h_his=False,
    bound_chainid_dict=None,
    all_models=False,
    update_coords=False,
    num_cpu=1,
):

    append_path(get_dir_path(dir_str=core_str, dir_path=core_dir))

    pdb_id_lst = type_lst(pdb_ids)

    pdb_dict = build_pdb_dict(pdb_id_lst)

    prot_resid_lst = res_to_lst(prot_resids)
    lig_resid_lst = res_to_lst(lig_resids)

    df = pd.DataFrame()

    if num_cpu == 1:
        for pdb_code in tqdm(
            list(pdb_dict.keys()),
            desc="Isolating PDB chains",
            position=0,
            leave=True,
        ):
            df = pd.concat(
                [
                    df,
                    isolate_chains(
                        pdb_code,
                        pdb_dict,
                        core_dir=core_dir,
                        rcsb_dir=rcsb_dir,
                        sifts_dir=sifts_dir,
                        renum_dir=renum_dir,
                        max_lig_dist=max_lig_dist,
                        min_lig_cont=min_lig_cont,
                        prot_resid_lst=prot_resid_lst,
                        lig_resid_lst=lig_resid_lst,
                        add_h=add_h,
                        add_h_his=add_h_his,
                        bound_chainid_dict=bound_chainid_dict,
                        all_models=all_models,
                        update_coords=update_coords,
                    ),
                ],
                sort=False,
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    isolate_chains,
                    pdb_code,
                    pdb_dict,
                    core_dir=core_dir,
                    rcsb_dir=rcsb_dir,
                    sifts_dir=sifts_dir,
                    renum_dir=renum_dir,
                    max_lig_dist=max_lig_dist,
                    min_lig_cont=min_lig_cont,
                    prot_resid_lst=prot_resid_lst,
                    lig_resid_lst=lig_resid_lst,
                    add_h=add_h,
                    add_h_his=add_h_his,
                    bound_chainid_dict=bound_chainid_dict,
                    all_models=all_models,
                    update_coords=update_coords,
                )
                for pdb_code in list(pdb_dict.keys())
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Isolating PDB chains",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                df = pd.concat([df, job.result()], sort=False)

    df = df.reset_index(drop=True)

    return df


def prep_coord(
    pdb_id_lst,
    renum_script_path,
    coord_table_path,
    core_dir=None,
    rcsb_dir=None,
    sifts_dir=None,
    renum_dir=None,
    sifts_json_path=None,
    max_lig_dist=4,
    min_lig_cont=5,
    prot_resids=None,
    lig_resids=None,
    add_h=False,
    add_h_his=False,
    bound_chainid_dict=None,
    all_models=False,
    data=None,
    update_coords=False,
    num_cpu=1,
):
    pdb_code_lst = build_pdb_code_lst(pdb_id_lst)

    run_pdb_renum(
        pdb_code_lst,
        renum_script_path,
        rcsb_dir=rcsb_dir,
        sifts_dir=sifts_dir,
        renum_dir=renum_dir,
        num_cpu=num_cpu,
    )

    delete_path("log_corrected.txt")
    delete_path("log_translator.txt")

    sifts_path_lst = [
        get_sifts_path(pdb_code, dir_path=sifts_dir) for pdb_code in pdb_code_lst
    ]

    build_sifts_map(sifts_path_lst, sifts_json_path, num_cpu=num_cpu)

    df = run_pdb_chain(
        pdb_id_lst,
        core_dir=core_dir,
        rcsb_dir=rcsb_dir,
        sifts_dir=sifts_dir,
        renum_dir=renum_dir,
        max_lig_dist=max_lig_dist,
        min_lig_cont=min_lig_cont,
        prot_resids=prot_resids,
        lig_resids=lig_resids,
        add_h=add_h,
        add_h_his=add_h_his,
        bound_chainid_dict=bound_chainid_dict,
        all_models=all_models,
        update_coords=update_coords,
        num_cpu=num_cpu,
    )

    if data is not None:
        df_col_lst = list(data.columns)

        for col in [pdb_id_col, bound_lig_col, bound_prot_chainid_col]:
            if col in df_col_lst:
                del data[col]

        df = merge_tables(df, data)

    save_table(coord_table_path, df)

    print("Prepared coordinate files!")
