# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""
from datetime import datetime

from .scripts import *
from .constants import *


def update_prep(data_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)
    sifts_json_path = get_file_path(sifts_json_file, dir_path=data_path)
    edia_json_path = get_file_path(edia_json_file, dir_path=data_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=data_path)

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
        delete_path(edia_json_path)
        delete_path(dih_json_path)

    if len(df) == 0:
        if past_df is not None:
            save_table(entry_table_path, past_df)
    else:
        pdb_id_lst = lst_col(df, pdb_id_col, unique=True)

        past_sifts_dict = load_json(sifts_json_path)

        prep_coord(
            pdb_id_lst=pdb_id_lst,
            renum_script_path=f"{get_dir_name(__file__)}/PDBrenum/PDBrenum.py",
            coord_table_path=entry_table_path,
            core_dir=data_path,
            rcsb_dir=data_path,
            sifts_dir=data_path,
            renum_dir=data_path,
            sifts_json_path=sifts_json_path,
            data=df,
            lig_resids="3-164",
            update_coords=False,
            num_cpu=num_cpu,
        )

        sifts_dict = load_json(sifts_json_path)
        if past_sifts_dict is not None:
            sifts_dict = merge_dicts([sifts_dict, past_sifts_dict])
        save_json(sifts_json_path, sifts_dict)

        df = load_table(entry_table_path)

        pdb_code_lst = lst_col(df, pdb_code_col, unique=True)
        sifts_dict = load_json(sifts_json_path)
        edia_dict = prep_edia(
            pdb_codes=pdb_code_lst,
            edia_dir=data_path,
            sifts_dict=sifts_dict,
            num_cpu=num_cpu,
        )

        past_edia_dict = load_json(edia_json_path)
        if past_edia_dict is not None:
            edia_dict = merge_dicts([edia_dict, past_edia_dict])
        save_json(edia_json_path, edia_dict)

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
                save_table(entry_table_path, df)


def update_annot(data_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)
    pdbaa_fasta_path = get_file_path(pdbaa_fasta_file, dir_path=data_path)

    df = load_table(entry_table_path)

    if past_df is not None:
        if (
            len(
                [
                    x
                    for x in [mut_class_col, pharm_class_col, prot_class_col, cf_col]
                    if x not in list(past_df.columns)
                ]
            )
            == 0
        ):
            df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))

    if len(df) == 0:
        if past_df is not None:
            save_table(entry_table_path, past_df)
    else:
        df = annot_mut(
            df=df,
            uniprot_accs=uniprot_acc_lst,
            resids="1-166",
            seq_dir=data_path,
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
            site_dict=pharm_site_dict,
            match_dict=pharm_match_dict,
            num_cpu=num_cpu,
        )

        for index in list(df.index.values):
            pharm_class = df.at[index, pharm_lig_site_col]
            if pharm_class not in [sp2_name, sp12_name, none_pharm_name]:
                pharm_class = other_pharm_name
            df.at[index, pharm_class_col] = pharm_class

        df[gene_class_col] = df[prot_col].map(gene_class_dict)
        df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

        df = annot_prot(
            df=df,
            pdbaa_fasta_path=pdbaa_fasta_path,
            pfam_dict=prot_pfam_dict,
        )

        for index in list(df.index.values):
            prot_class = df.at[index, bound_prot_pfam_col]
            for prot_name, pfam_lst in prot_class_dict.items():
                if prot_class in pfam_lst:
                    prot_class = prot_name
            if prot_class == gef_name:
                nuc_class = df.at[index, nuc_class_col]
                prot_class += "."
                if nuc_class == gtp_name:
                    prot_class += gef_rem_name
                else:
                    prot_class += gef_cdc_name
            df.at[index, prot_class_col] = prot_class

        if (
            past_df is not None
            and len([x for x in class_col_lst if x not in list(past_df.columns)]) == 0
        ):
            df[complete_col] = str(False)
            df = pd.concat([df, past_df], sort=False)
            save_table(entry_table_path, df)

        rcsb_path_lst = lst_col(df, rcsb_path_col, unique=True)
        df = annot_cf(
            coord_paths=rcsb_path_lst,
            data=df,
            num_cpu=num_cpu,
        )

        save_table(entry_table_path, df)


def update_interf(data_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)
    interf_json_path = get_file_path(interf_json_file, dir_path=data_path)
    interf_table_path = get_file_path(interf_table_file, dir_path=data_path)

    df = load_table(entry_table_path)

    delete_paths = True
    if past_df is not None and interf_class_col in list(past_df.columns):
        df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))
        delete_paths = False

    if delete_paths:
        delete_path(interf_json_path)
        delete_path(interf_table_path)

    if len(df) == 0:
        if past_df is not None:
            save_table(entry_table_path, past_df)
    else:
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
            interf_dir=data_path,
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
            search_coord_path=get_renum_path(sup_pdb_code, dir_path=data_path),
            search_chainid=sup_chainid,
            search_interf="4",
        )

        past_interf_df = load_table(interf_table_path)
        if past_interf_df is not None and interf_class_col in list(past_df.columns):
            interf_df = pd.concat([interf_df, past_interf_df], sort=False)
        save_table(interf_table_path, interf_df)

        interf_pdb_id_lst = lst_col(interf_df, pdb_id_col)
        df[interf_class_col] = (
            df[pdb_id_col]
            .map(make_dict(interf_pdb_id_lst, len(interf_pdb_id_lst) * [dimer_name]))
            .fillna(none_dimer_name)
        )

        if past_df is not None:
            if interf_class_col in list(past_df.columns):
                df[complete_col] = str(False)
                df = pd.concat([df, past_df], sort=False)

        save_table(entry_table_path, df)


def update_pocket(data_path=None, past_df=None, num_cpu=1):

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)
    pocket_json_path = get_file_path(pocket_json_file, dir_path=data_path)
    pocket_table_path = get_file_path(pocket_table_file, dir_path=data_path)

    df = load_table(entry_table_path)

    delete_paths = True
    if past_df is not None:
        if pocket_class_col in list(past_df.columns):
            df = mask_unequal(df, pdb_id_col, lst_col(past_df, pdb_id_col))
            delete_paths = False

    if delete_paths:
        delete_path(pocket_json_path)
        delete_path(pocket_table_path)

    if len(df) == 0:
        if past_df is not None:
            save_table(entry_table_path, past_df)
    else:
        coord_path_lst = [x.replace(".cif", ".pdb") for x in lst_col(df, core_path_col)]

        pharm_dict = dict()
        chainid_dict = dict()
        for coord_path in coord_path_lst:
            chainid_dict[coord_path] = list()
        for index in list(df.index.values):
            coord_path = df.at[index, core_path_col].replace(".cif", ".pdb")
            chainid = df.at[index, chainid_col]
            chainid_dict[coord_path].append(chainid)
            if df.at[index, pharm_lig_site_col] in [sp2_name, sp12_name]:
                pharm_dict[coord_path] = df.at[index, pharm_lig_col]

        pocket_dict = prep_pocket(
            coord_paths=coord_path_lst,
            pocket_dir=data_path,
            pharm_dict=pharm_dict,
            chainid_dict=chainid_dict,
            num_cpu=num_cpu,
        )

        past_pocket_dict = load_json(pocket_json_path)
        if past_pocket_dict is not None:
            pocket_dict = merge_dicts([pocket_dict, past_pocket_dict])
        save_json(pocket_json_path, pocket_dict)

        if past_df is not None:
            if pocket_class_col in list(past_df.columns):
                df[complete_col] = str(False)
                df = pd.concat([df, past_df], sort=False)

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

            search_cont_lst = lst_col(
                mask_equal(bound_df, pharm_class_col, pharm_name), pharm_lig_cont_col
            )
            search_cont_lst = [str_to_lst(x) for x in search_cont_lst]

            unbound_df = build_pocket_table(
                df=mask_equal(df, pharm_class_col, none_pharm_name),
                pocket_dict=pocket_dict,
                search_cont_lst=search_cont_lst,
                search_max_dist=0.4,
                use_simpson=True,
            )

            unbound_df = mask_equal(unbound_df, pocket_status_col, pocket_unbound_name)

            bound_df[pocket_site_col] = pharm_name
            unbound_df[pocket_site_col] = pharm_name

            pocket_df = pd.concat([pocket_df, bound_df], sort=False)
            pocket_df = pd.concat([pocket_df, unbound_df], sort=False)

        other_df = mask_unequal(
            temp_df, pocket_id_col, lst_col(pocket_df, pocket_id_col)
        )
        other_df = mask_equal(other_df, pocket_status_col, pocket_unbound_name)
        other_df = mask_equal(other_df, pharm_class_col, none_pharm_name)

        other_df[pocket_site_col] = other_pharm_name

        pocket_df = pd.concat([pocket_df, other_df], sort=False)

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
                if other_pharm_name in pocket_site_lst:
                    pocket_site_lst.remove(other_pharm_name)

            if len(pocket_site_lst) == 0:
                pharm_class = df.at[index, pharm_class_col]
                if pharm_class in [sp2_name, sp12_name, other_pharm_name]:
                    pocket_site_lst = [pharm_class]
                    pocket_status_lst = [pocket_bound_name]
                else:
                    pocket_site_lst = [other_pharm_name]
                    pocket_status_lst = [pocket_unbound_name]

            pocket_site = lst_to_str(pocket_site_lst, join_txt="|")
            pocket_status = lst_to_str(pocket_status_lst)

            df.at[index, pocket_site_col] = pocket_site
            df.at[index, pocket_status_col] = pocket_status

            pocket_class = pocket_site

            pocket_class += "-"
            pocket_class += pocket_status

            df.at[index, pocket_class_col] = pocket_class

        save_table(entry_table_path, df)


def update_rascore(data_path=None, pdbaa_fasta_path=None, num_cpu=1):

    if data_path is None:
        data_path = f"{os.getcwd()}/{rascore_str}_{data_str}"

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)

    if pdbaa_fasta_path is None:
        curr_date = datetime.today().strftime("%Y-%m-%d")
        pdbaa_fasta_path = get_file_path(
            f"{curr_date}_{pdbaa_fasta_file}",
            dir_str=pdbaa_str,
            dir_path=data_path,
            pre_str=False,
        )
    else:
        delete_path(entry_table_path)

    past_df = load_table(entry_table_path)

    if past_df is None:
        print("No rascore database found! Building database from scratch!")
    else:
        if complete_col in list(past_df.columns):
            past_df = mask_equal(past_df, complete_col, str(True))
        else:
            past_df[complete_col] = str(True)

    print("Downloading updated pdbaa file.")
    download_unzip(url=pdbaa_url, path=pdbaa_fasta_path)

    search_pdbaa(
        pdbaa_fasta_path=pdbaa_fasta_path,
        search_lst=swiss_id_lst,
        entry_table_path=entry_table_path,
        min_length=25,
    )

    if past_df is not None:
        df = load_table(entry_table_path)
        if len(df) < len(past_df):
            print("Downgrading rascore database to older version.")
            past_df = None

    update_prep(data_path=data_path, past_df=past_df, num_cpu=num_cpu)
    update_annot(data_path=data_path, past_df=past_df, num_cpu=num_cpu)
    update_interf(data_path=data_path, past_df=past_df, num_cpu=num_cpu)
    update_pocket(data_path=data_path, past_df=past_df, num_cpu=num_cpu)

    df = load_table(entry_table_path)
    if complete_col in list(df.columns):
        del df[complete_col]
    save_table(entry_table_path, df)

    print("Rascore update complete!")