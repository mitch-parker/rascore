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
import pandas as pd
from tqdm import tqdm
from datetime import datetime

from ..pipelines.classify_rascore import classify_rascore

from ..constants.conf import sw1_name, sw2_name, y32_name, y71_name
from ..constants.nuc import nuc_class_dict, gtp_name
from ..constants.gene import uniprot_acc_lst, swiss_id_lst, gene_class_dict
from ..constants.mut import mut_class_lst, other_mut_name
from ..constants.pharm import (
    pharm_match_dict,
    pharm_site_dict,
    sp2_name,
    sp12_name,
    mult_pharm_name,
    none_pharm_name,
    other_pharm_name,
)
from ..constants.dimer import dimer_name, none_dimer_name
from ..constants.prot import (
    gef_name,
    prot_pfam_dict,
    prot_class_dict,
    rem_name,
    cdc_name,
    nano_name,
    mult_prot_name,
)
from ..constants.pml import sup_pdb_code, sup_chainid


from ..scripts.search_pdbaa import search_pdbaa
from ..scripts.prep_coord import prep_coord
from ..scripts.prep_dih import prep_dih
from ..scripts.annot_mut import annot_mut
from ..scripts.annot_lig import annot_lig
from ..scripts.annot_prot import annot_prot
from ..scripts.annot_cf import annot_cf
from ..scripts.prep_interf import prep_interf
from ..scripts.build_interf_table import build_interf_table
from ..scripts.prep_pocket import prep_pocket, pocket_bound_name, pocket_unbound_name
from ..scripts.build_pocket_table import build_pocket_table

from ..functions.table import (
    mask_equal,
    make_dict,
    merge_dicts,
    mask_unequal,
    lst_col,
    title_str,
)
from ..functions.path import (
    get_renum_path,
    get_dir_path,
    copy_path,
    load_json,
    path_exists,
    save_table,
    save_json,
    get_file_path,
    load_table,
    get_neighbor_path,
    delete_path,
    get_file_name,
    pipelines_str,
    data_str,
    classify_str,
    rascore_str,
    build_str,
    pdbaa_str,
)
from ..functions.lst import (
    lst_to_str,
    str_to_lst,
)
from ..functions.file import (
    entry_table_file,
    sifts_json_file,
    dih_json_file,
    interf_json_file,
    interf_table_file,
    pocket_json_file,
    pocket_table_file,
    result_table_file,
    pdbaa_fasta_file,
)
from ..functions.url import pdbaa_url
from ..functions.download import download_unzip
from ..functions.col import (
    pdb_id_col,
    chainid_col,
    core_path_col,
    bound_lig_col,
    bound_prot_chainid_col,
    mut_class_col,
    pharm_class_col,
    prot_class_col,
    cf_col,
    mut_status_col,
    pharm_lig_site_col,
    pharm_lig_col,
    pharm_lig_match_col,
    match_class_col,
    gene_class_col,
    prot_class_col,
    prot_col,
    bio_lig_col,
    nuc_class_col,
    rcsb_path_col,
    bound_prot_pfam_col,
    complete_col,
    interf_class_col,
    method_col,
    method_col,
    renum_path_col,
    pocket_class_col,
    pocket_type_col,
    pocket_status_col,
    pocket_site_col,
    pocket_id_col,
    pocket_cont_col,
)


def update_prep(out_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)
    sifts_json_path = get_file_path(sifts_json_file, dir_path=out_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=out_path)

    df = load_table(entry_table_path)

    delete_paths = True
    if past_df is not None:
        if (
            len(
                [
                    x
                    for x in [bound_lig_col, bound_prot_chainid_col]
                    if x not in list(past_df.columns)
                ]
            )
            == 0
        ):
            df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))
            delete_paths = False

    if delete_paths:
        delete_path(sifts_json_path)
        delete_path(dih_json_path)
        past_df = None

    if len(df) > 0:
        pdb_id_lst = lst_col(df, pdb_id_col, unique=True)

        past_sifts_dict = load_json(sifts_json_path)

        prep_coord(
            pdb_id_lst=pdb_id_lst,
            renum_script_path=get_neighbor_path(__file__, pipelines_str, "PDBrenum")
            + "/PDBrenum.py",
            coord_table_path=entry_table_path,
            core_dir=out_path,
            rcsb_dir=out_path,
            sifts_dir=out_path,
            renum_dir=out_path,
            sifts_json_path=sifts_json_path,
            data=df,
            update_coords=False,
            num_cpu=num_cpu,
        )

        sifts_dict = load_json(sifts_json_path)
        if past_sifts_dict is not None:
            sifts_dict = merge_dicts([sifts_dict, past_sifts_dict])
        save_json(sifts_json_path, sifts_dict)
        
        try:
            df = load_table(entry_table_path)

            core_path_lst = lst_col(df, core_path_col, unique=True)
            dih_dict = prep_dih(coord_paths=core_path_lst, num_cpu=num_cpu)

            past_dih_dict = load_json(dih_json_path)
            if past_dih_dict is not None:
                dih_dict = merge_dicts([dih_dict, past_dih_dict])
            save_json(dih_json_path, dih_dict)

            if past_df is not None:
                if (
                    len(
                        [
                            x
                            for x in [bound_lig_col, bound_prot_chainid_col]
                            if x not in list(past_df.columns)
                        ]
                    )
                    == 0
                ):
                    df = pd.concat([df, past_df], sort=False)
                    df = df.reset_index(drop=True)
                    save_table(entry_table_path, df)
        except:
            if past_df is not None:
                save_table(entry_table_path, past_df)



def update_annot(pdbaa_fasta_path, out_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)

    df = load_table(entry_table_path)

    if past_df is not None:
        if (
            len(
                [
                    x
                    for x in [
                        gene_class_col,
                        nuc_class_col,
                        mut_class_col,
                        pharm_class_col,
                        match_class_col,
                        prot_class_col,
                        cf_col,
                    ]
                    if x not in list(past_df.columns)
                ]
            )
            == 0
        ):
            df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))

    if len(df) > 0:
        df = annot_mut(
            df=df,
            uniprot_accs=uniprot_acc_lst,
            resids="1-166",
            seq_dir=out_path,
            num_cpu=num_cpu,
        )

        for index in list(df.index.values):
            df.at[index, mut_class_col] = lst_to_str(
                [
                    x
                    for x in mut_class_lst
                    if x in str_to_lst(df.at[index, mut_status_col])
                ],
                empty=other_mut_name,
            )

        df = annot_lig(
            df=df,
            lig_dir=out_path,
            site_dict=pharm_site_dict,
            match_dict=pharm_match_dict,
            num_cpu=num_cpu,
        )

        for index in list(df.index.values):
            pharm_class = other_pharm_name
            pharm_site = df.at[index, pharm_lig_site_col]
            if "|" in pharm_site:
                pharm_class = mult_pharm_name
            elif pharm_site in [sp2_name, sp12_name, none_pharm_name]:
                pharm_class = pharm_site
            df.at[index, pharm_class_col] = pharm_class

        for index in list(df.index.values):
            pharm_match = df.at[index, pharm_lig_match_col]
            if "|" in pharm_match:
                match_class = mult_pharm_name
            else:
                match_class = pharm_match
            df.at[index, match_class_col] = match_class

        df[gene_class_col] = df[prot_col].map(gene_class_dict)
        df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

        df = annot_prot(
            df=df,
            pdbaa_fasta_path=pdbaa_fasta_path,
            pfam_dict=prot_pfam_dict,
        )

        for index in list(df.index.values):
            prot_pfam = df.at[index, bound_prot_pfam_col]
            if "|" in prot_pfam:
                if nano_name in prot_pfam:
                    prot_class = nano_name
                else:
                    prot_class = mult_prot_name
            else:
                for prot_name, pfam_lst in prot_class_dict.items():
                    if prot_pfam in pfam_lst:
                        prot_class = prot_name
                if prot_class == gef_name:
                    nuc_class = df.at[index, nuc_class_col]
                    prot_class += "."
                    if nuc_class == gtp_name:
                        prot_class += rem_name
                    else:
                        prot_class += cdc_name
            df.at[index, prot_class_col] = prot_class

        if (
            past_df is not None
            and len([x for x in [
                        gene_class_col,
                        nuc_class_col,
                        mut_class_col,
                        pharm_class_col,
                        match_class_col,
                        prot_class_col,
                        cf_col,
                    ] if x not in list(past_df.columns)]) == 0
        ):
            df[complete_col] = str(False)
            df = pd.concat([df, past_df], sort=False)
            df = df.reset_index(drop=True)
            save_table(entry_table_path, df)

        rcsb_path_lst = lst_col(df, rcsb_path_col, unique=True)
        df = annot_cf(
            coord_paths=rcsb_path_lst,
            data=df,
            num_cpu=num_cpu,
        )

        save_table(entry_table_path, df)


def update_interf(out_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)
    interf_json_path = get_file_path(interf_json_file, dir_path=out_path)
    interf_table_path = get_file_path(interf_table_file, dir_path=out_path)

    df = load_table(entry_table_path)

    delete_paths = True
    if past_df is not None and interf_class_col in list(past_df.columns) and path_exists(interf_json_path) and path_exists(interf_table_path):
        df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))
        delete_paths = False

    if delete_paths:
        delete_path(interf_json_path)
        delete_path(interf_table_path)
        past_df = None

    if len(df) > 0:
        xray_df = mask_equal(df, method_col, "XRAY")

        coord_path_lst = lst_col(xray_df, renum_path_col, unique=True)

        chainid_dict = dict()
        for coord_path in coord_path_lst:
            chainid_dict[coord_path] = list()
        for index in list(xray_df.index.values):
            coord_path = xray_df.at[index, renum_path_col]
            chainid = xray_df.at[index, chainid_col]
            chainid_dict[coord_path].append(chainid)

        interf_dict = prep_interf(
            coord_paths=coord_path_lst,
            interf_dir=out_path,
            chainid_dict=chainid_dict,
            num_cpu=num_cpu,
        )

        past_interf_dict = load_json(interf_json_path)
        if past_interf_dict is not None:
            interf_dict = merge_dicts([interf_dict, past_interf_dict])
        save_json(interf_json_path, interf_dict)

        interf_df = build_interf_table(
            df=xray_df,
            interf_dict=interf_dict,
            search_coord_path=get_renum_path(sup_pdb_code, dir_path=out_path),
            search_chainid=sup_chainid,
            search_interf="4",
        )

        past_interf_df = load_table(interf_table_path)
        if past_interf_df is not None and interf_class_col in list(past_df.columns):
            interf_df = pd.concat([interf_df, past_interf_df], sort=False)
            interf_df = interf_df.reset_index(drop=True)
        save_table(interf_table_path, interf_df)

        interf_pdb_id_lst = lst_col(interf_df, pdb_id_col)
        df[interf_class_col] = (
            df[pdb_id_col]
            .map(make_dict(interf_pdb_id_lst, len(interf_pdb_id_lst) * [dimer_name]))
            .fillna(none_dimer_name)
        )

        if past_df is not None and interf_class_col in list(past_df.columns):
            df[complete_col] = str(False)
            df = pd.concat([df, past_df], sort=False)
            df = df.reset_index(drop=True)

        save_table(entry_table_path, df)


def update_pocket(out_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)
    pocket_json_path = get_file_path(pocket_json_file, dir_path=out_path)
    pocket_table_path = get_file_path(pocket_table_file, dir_path=out_path)

    df = load_table(entry_table_path)

    delete_paths = True
    if past_df is not None and pocket_class_col in list(past_df.columns) and path_exists(pocket_json_path) and path_exists(pocket_table_path):
        df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))
        delete_paths = False

    if delete_paths:
        delete_path(pocket_json_path)
        delete_path(pocket_table_path)
        past_df = None

    if len(df) > 0:
        coord_path_lst = [x.replace(".cif", ".pdb") for x in lst_col(df, core_path_col)]

        chainid_dict = dict()
        for coord_path in coord_path_lst:
            chainid_dict[coord_path] = list()
        for index in list(df.index.values):
            coord_path = df.at[index, core_path_col].replace(".cif", ".pdb")
            chainid = df.at[index, chainid_col]
            chainid_dict[coord_path].append(chainid)

        pocket_dict = prep_pocket(
            coord_paths=coord_path_lst,
            pocket_dir=out_path,
            chainid_dict=chainid_dict,
            num_cpu=num_cpu,
        )

        past_pocket_dict = load_json(pocket_json_path)
        if past_pocket_dict is not None:
            pocket_dict = merge_dicts([pocket_dict, past_pocket_dict])
        save_json(pocket_json_path, pocket_dict)

        if past_df is not None and pocket_class_col in list(past_df.columns):
            df[complete_col] = str(False)
            df = pd.concat([df, past_df], sort=False)
            df = df.reset_index(drop=True)

        temp_df = build_pocket_table(
            df=df,
            pocket_dict=pocket_dict,
        )
        temp_df = mask_equal(
            temp_df, pocket_type_col, [title_str(pharm_lig_col), pocket_unbound_name]
        )

        pocket_df = pd.DataFrame()
        for pharm_name in [sp2_name, sp12_name]:
            bound_df = mask_equal(temp_df, pharm_class_col, pharm_name)
            bound_df = mask_equal(bound_df, pocket_status_col, pocket_bound_name)

            bound_df[pocket_type_col] = pharm_name
            bound_df[pocket_site_col] = pharm_name
            pocket_df = pd.concat([pocket_df, bound_df], sort=False)
            pocket_df = pocket_df.reset_index(drop=True)

            search_cont_lst = lst_col(bound_df, pocket_cont_col)
            search_cont_lst = [str_to_lst(x) for x in search_cont_lst]

            unbound_df = build_pocket_table(
                df=mask_equal(df, pharm_class_col, none_pharm_name),
                pocket_dict=pocket_dict,
                search_cont_lst=search_cont_lst,
                search_max_dist=0.4,
                use_simpson=True,
            )

            if len(unbound_df) > 0:
                unbound_df = mask_equal(
                    unbound_df, pocket_status_col, pocket_unbound_name
                )
                unbound_df[pocket_site_col] = pharm_name
                pocket_df = pd.concat([pocket_df, unbound_df], sort=False)
                pocket_df = pocket_df.reset_index(drop=True)

        other_df = mask_unequal(
            temp_df, pocket_id_col, lst_col(pocket_df, pocket_id_col)
        )
        other_df = mask_equal(other_df, pocket_status_col, pocket_unbound_name)
        other_df = mask_equal(other_df, pharm_class_col, none_pharm_name)

        other_df[pocket_site_col] = other_pharm_name

        pocket_df = pd.concat([pocket_df, other_df], sort=False)
        pocket_df = pocket_df.reset_index(drop=True)

        save_table(pocket_table_path, pocket_df)

        pocket_df = load_table(pocket_table_path)

        annot_df = mask_equal(pocket_df, pocket_site_col, [sp2_name, sp12_name])

        for index in tqdm(
            list(df.index.values),
            desc="Annotating pockets",
            position=0,
            leave=True,
        ):
            pdb_df = mask_equal(annot_df, pdb_id_col, df.at[index, pdb_id_col])

            pocket_site_lst = lst_col(pdb_df, pocket_site_col, unique=True)
            pocket_status_lst = lst_col(pdb_df, pocket_status_col, unique=True)

            if len(pocket_site_lst) > 1:
                pocket_site_lst = [mult_pharm_name]

            if len(pocket_site_lst) == 0:
                pharm_class = df.at[index, pharm_class_col]
                if pharm_class == none_pharm_name:
                    pocket_site_lst = [other_pharm_name]
                    pocket_status_lst = [pocket_unbound_name]
                else:
                    pocket_site_lst = [pharm_class]
                    pocket_status_lst = [pocket_bound_name]

            pocket_site = lst_to_str(pocket_site_lst, join_txt="|")
            pocket_status = lst_to_str(pocket_status_lst)

            df.at[index, pocket_site_col] = pocket_site
            df.at[index, pocket_status_col] = pocket_status

            pocket_class = pocket_site

            pocket_class += "-"
            pocket_class += pocket_status

            df.at[index, pocket_class_col] = pocket_class

        save_table(entry_table_path, df)


def update_classify(out_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=out_path)

    df = load_table(entry_table_path)

    if past_df is not None:
        if (
            len(
                [
                    x
                    for x in [
                        sw1_name,
                        sw2_name,
                        y32_name,
                        y71_name
                    ]
                    if x not in list(past_df.columns)
                ]
            )
            == 0
        ):
            df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))

    dih_dict = load_json(dih_json_path)

    classify_path = f"{out_path}/{classify_str}"

    if len(df) > 0:
        classify_rascore(
            df,
            out_path=classify_path,
            dih_dict=dih_dict,
            num_cpu=num_cpu,
        )

        result_df = load_table(get_file_path(result_table_file, dir_path=classify_path))

        for col in [sw1_name, sw2_name, y32_name, y71_name]:
            df[col] = df[pdb_id_col].map(
                make_dict(lst_col(result_df, pdb_id_col), lst_col(result_df, col))
            )

        if past_df is not None:
            if (
                len(
                    [
                        x
                        for x in [
                            sw1_name,
                            sw2_name,
                            y32_name,
                            y71_name
                        ]
                        if x not in list(past_df.columns)
                    ]
                )
                == 0
            ):
                df[complete_col] = str(False)
                df = pd.concat([df, past_df], sort=False)
                df = df.reset_index(drop=True)

        save_table(entry_table_path, df)


def build_rascore(out_path=None, pdbaa_fasta_path=None, num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{build_str}"

    entry_table_path = get_file_path(entry_table_file, dir_path=out_path)
    pocket_table_path = get_file_path(pocket_table_file, dir_path=out_path)
    interf_table_path = get_file_path(interf_table_file, dir_path=out_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=out_path)

    if pdbaa_fasta_path is None:
        curr_date = datetime.today().strftime("%Y-%m-%d")
        pdbaa_fasta_path = get_file_path(
            f"{curr_date}_{pdbaa_fasta_file}",
            dir_str=pdbaa_str,
            dir_path=out_path,
            pre_str=False,
        )

    past_df = load_table(entry_table_path)

    if past_df is None:
        print("No rascore database found! Building rascore database from scratch!")
    else:
        if complete_col in list(past_df.columns):
            past_df = mask_equal(past_df, complete_col, str(True))
        else:
            past_df[complete_col] = str(True)

    print("Downloading updated pdbaa file.")

    try:
        download_unzip(url=pdbaa_url, path=pdbaa_fasta_path)
    except:
        pdbaa_fasta_path_lst = [
            x
            for x in os.listdir(get_dir_path(dir_str=pdbaa_str, dir_path=out_path))
            if pdbaa_fasta_file in x
        ]
        if len(pdbaa_fasta_path_lst) > 0:
            pdbaa_fasta_path = sorted(
                pdbaa_fasta_path_lst,
                reverse=True,
            )[0]

            last_date = pdbaa_fasta_path.split(f"_{pdbaa_fasta_file}")[0]
            last_date = last_date.split(f"{pdbaa_str}/")[1]
            print(f"Updated pdbaa file unavailable. Using latest version ({last_date})")

    try:
        search_pdbaa(
            pdbaa_fasta_path=pdbaa_fasta_path,
            search_lst=swiss_id_lst,
            entry_table_path=entry_table_path,
            min_length=25,
        )
    except:
        print("Error reading pdbaa file.")
        copy_path(
            get_file_path(
                entry_table_file,
                dir_path=get_neighbor_path(__file__, pipelines_str, data_str),
            ),
            entry_table_path,
        )

    if past_df is not None:
        df = load_table(entry_table_path)
        if len(df) < len(past_df):
            print("Downgrading rascore database to older version.")
            past_df = None

    update_prep(out_path=out_path, past_df=past_df, num_cpu=num_cpu)
    update_annot(
        pdbaa_fasta_path,
        out_path=out_path,
        past_df=past_df,
        num_cpu=num_cpu,
    )
    update_interf(out_path=out_path, past_df=past_df, num_cpu=num_cpu)
    update_pocket(out_path=out_path, past_df=past_df, num_cpu=num_cpu)
    update_classify(out_path=out_path, past_df=past_df, num_cpu=num_cpu)

    df = load_table(entry_table_path)
    if complete_col in list(df.columns):
        del df[complete_col]

    save_table(entry_table_path, df)

    ref_df = df.set_index(pdb_id_col)

    pocket_df = load_table(pocket_table_path)
    interf_df = load_table(interf_table_path)

    for col in list(ref_df.columns):
        if col not in list(pocket_df.columns):
            for index in list(pocket_df.index.values):
                pocket_df.at[index, col] = ref_df.at[pocket_df.at[index, pdb_id_col], col]
        if col not in list(interf_df.columns):
            for index in list(interf_df.index.values):
                interf_df.at[index, col] = ref_df.at[interf_df.at[index, pdb_id_col], col]

    save_table(pocket_table_path, pocket_df)
    save_table(interf_table_path, interf_df)

    for file_path in [
        entry_table_path,
        pocket_table_path,
        interf_table_path,
        dih_json_path,
    ]:
        print(get_file_path(
                get_file_name(file_path),
                dir_path=get_neighbor_path(__file__, pipelines_str, data_str),
            ))
        copy_path(
            file_path,
            get_file_path(
                get_file_name(file_path),
                dir_path=get_neighbor_path(__file__, pipelines_str, data_str),
            ),
        )

    print("Rascore update complete!")