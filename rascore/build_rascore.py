# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from scripts import *
from constants import *


def build_rascore(pdbaa_fasta_path=None, data_path=None, num_cpu=1):

    if data_path is None:
        data_path = get_dir_path(dir_str=data_str, dir_path=get_dir_name(__file__))

    if pdbaa_fasta_path is None:
        pdbaa_fasta_path = get_file_path(pdbaa_fasta_file, dir_path=data_path)
        download_unzip(url=pdbaa_url, path=pdbaa_fasta_path)

    entry_table_path = get_file_path(entry_table_file, dir_path=data_path)

    search_pdbaa(
        pdbaa_fasta_path=pdbaa_fasta_path,
        search_lst=swiss_id_lst,
        entry_table_path=entry_table_path,
        min_length=25,
    )
    df = load_table(entry_table_path)

    pdb_id_lst = lst_col(df, pdb_id_col, unique=True)

    sifts_json_path = get_file_path(sifts_json_file, dir_path=data_path)
    prep_coord(
        pdb_id_lst=pdb_id_lst,
        renum_script_path=f"{get_dir_name(__file__)}/scripts/PDBrenum/PDBrenum.py",
        coord_table_path=entry_table_path,
        sifts_json_path=sifts_json_path,
        data=df,
        num_cpu=num_cpu,
    )
    df = load_table(entry_table_path)

    annot_mut(
        df=df,
        uniprot_accs=uniprot_acc_lst,
        resids="1-166",
        mut_table_path=entry_table_path,
        num_cpu=num_cpu,
    )
    df = load_table(entry_table_path)

    for index in list(df.index.values):
        df.at[index, mut_class_col] = lst_to_str(
            [x for x in mut_class_lst if x in str_to_lst(df.at[index, mut_status_col])],
            empty=other_mut_name,
        )
    save_table(entry_table_path, df)

    annot_lig(
        df=df,
        site_dict=pharm_site_dict,
        match_dict=pharm_match_dict,
        lig_table_path=entry_table_path,
        num_cpu=num_cpu,
    )
    df = load_table(entry_table_path)

    for index in list(df.index.values):
        pharm_status = df.at[index, pharm_lig_site_col]
        if pharm_status == sp2_name or pharm_status == sp12_name:
            pharm_status += "."
            pharm_status += df.at[index, pharm_lig_match_col]
        df.at[index, pharm_class_col] = pharm_status
    save_table(entry_table_path, df)

    annot_prot(
        df=df,
        pdbaa_fasta_path=pdbaa_fasta_path,
        pfam_dict=prot_pfam_dict,
        prot_table_path=entry_table_path,
    )
    df = load_table(entry_table_path)

    for index in list(df.index.values):
        prot_status = df.at[index, bound_prot_pfam_col]
        for prot_class, pfam_lst in prot_class_dict.items():
            if prot_status in pfam_lst:
                prot_status = prot_class
        if prot_status == gef_name:
            nuc_class = df.at[index, nuc_class_col]
            prot_class += "."
            if nuc_class == gtp_name:
                prot_class += gef_rem_name
            else:
                prot_class += gef_cdc_name
        df.at[index, prot_class_col] = prot_status
    save_table(entry_table_path, df)

    df[gene_class_col] = df[prot_col].map(gene_class_dict)
    df[nuc_class_col] = df[bio_lig_col].map(nuc_class_dict).fillna(gtp_name)

    rcsb_path_lst = lst_col(df, rcsb_path_col, unique=True)
    annot_cf(
        coord_paths=rcsb_path_lst,
        cf_table_path=entry_table_path,
        data=df,
        num_cpu=num_cpu,
    )
    df = load_table(entry_table_path)

    edia_json_path = get_file_path(edia_json_file, dir_path=data_path)
    pdb_code_lst = lst_col(df, pdb_code_col, unique=True)
    sifts_dict = load_json(sifts_json_path)
    prep_edia(
        pdb_codes=pdb_code_lst,
        edia_json_path=edia_json_path,
        sifts_dict=sifts_dict,
        num_cpu=num_cpu,
    )

    dih_json_path = get_file_path(dih_json_file, dir_path=data_path)
    core_path_lst = lst_col(df, core_path_col, unique=True)
    prep_dih(coord_paths=core_path_lst, dih_json_path=dih_json_path, num_cpu=num_cpu)
