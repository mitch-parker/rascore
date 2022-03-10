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

from ..functions.seq import load_record_dict, get_record_desc
from ..functions.lst import lst_to_str, res_to_lst, str_to_lst, lst_unique, sort_lst
from ..functions.coord import (
    load_coord,
    get_neighbors,
    get_resnum,
    get_chain_cont,
    get_reschainid,
    is_aa,
)
from ..functions.col import (
    core_path_col,
    modelid_col,
    chainid_col,
    pdb_code_col,
    bound_prot_chainid_col,
    bound_prot_col,
    bound_prot_site_col,
    bound_prot_pfam_col,
    bound_prot_swiss_id_col,
    bound_prot_cont_col
)
from ..functions.table import fix_val
from ..functions.path import save_table
from ..functions.pdbaa import get_pdbaa_prot, get_pdbaa_swiss_id


def annot_prot(
    df,
    pdbaa_fasta_path,
    pfam_dict=None,
    site_dict=None,
    other_pfam="Other",
    other_site="Unclassified",
    none_pfam="None",
    none_site="None",
    max_site_dist=4,
    prot_table_path=None,
):

    pdbaa_dict = load_record_dict(pdbaa_fasta_path)

    if site_dict is not None:
        for site_name, site_cont in site_dict.items():
            site_dict[site_name] = res_to_lst(site_cont)

    if pfam_dict is not None:
        pfam_dict_lst = list(pfam_dict.keys())

    for index in tqdm(
        list(df.index.values), desc="Annotating proteins", position=0, leave=True
    ):

        pdb_code = df.at[index, pdb_code_col]

        bound_prot_chainid = df.at[index, bound_prot_chainid_col]

        prot_lst = list()
        cont_lst = list()
        swiss_id_lst = list()
        site_lst = list()
        pfam_lst = list()

        if bound_prot_chainid != "None":
            coord_path = df.at[index, core_path_col]
            modelid = df.at[index, modelid_col]
            chainid = df.at[index, chainid_col]
            structure = load_coord(coord_path)
            neighbors = get_neighbors(
                structure[fix_val(modelid, return_int=True)][chainid]
            )
            for bound_chainid in str_to_lst(bound_prot_chainid):
                record = pdbaa_dict[f"{pdb_code.upper()}{bound_chainid}"]
                desc = get_record_desc(record)
                prot = get_pdbaa_prot(desc)
                swiss_id = get_pdbaa_swiss_id(desc)

                prot_lst.append(prot)
                swiss_id_lst.append(swiss_id)

                resid_cont_lst = lst_unique(
                    [
                        get_resnum(x)
                        for x in get_chain_cont(
                            neighbors,
                            structure[fix_val(modelid, return_int=True)][
                                bound_chainid
                            ],
                            level="R",
                            max_dist=max_site_dist,
                        )
                        if (get_resnum(x) < 50000)
                        and (get_reschainid(x) == chainid)
                        and is_aa(x)
                    ]
                )

                if len(resid_cont_lst) > 0:
                    resid_cont_str = bound_chainid
                    resid_cont_str += ':'
                    resid_cont_str += lst_to_str(resid_cont_lst)
                    cont_lst.append(resid_cont_str)

                if site_dict is not None:

                    max_cont = 0
                    prot_site = other_site
                    for site_name, site_cont_lst in site_dict.items():
                        site_cont = len([x for x in site_cont_lst if x in resid_cont_lst])
                        if site_cont > max_cont:
                            prot_site = site_name
                            max_cont = site_cont
                    site_lst.append(prot_site)

                if pfam_dict is not None:
                    pfam_status = other_pfam
                    if swiss_id in pfam_dict_lst:
                        pfam_status = pfam_dict[swiss_id]
                    elif swiss_id == "NA":
                        pfam_status = "Binder"
                    pfam_lst.append(pfam_status)

        df.at[index, bound_prot_col] = lst_to_str(sort_lst(lst_unique(prot_lst)))
        df.at[index, bound_prot_swiss_id_col] = lst_to_str(
            sort_lst(lst_unique(swiss_id_lst))
        )
        df.at[index, bound_prot_site_col] = lst_to_str(
            sort_lst(lst_unique(site_lst)), join_txt="|", empty=none_site
        )
        df.at[index, bound_prot_pfam_col] = lst_to_str(
            sort_lst(lst_unique(pfam_lst)), join_txt="|", empty=none_pfam
        )

        df.at[index, bound_prot_cont_col] = lst_to_str(cont_lst, join_txt="|")

    print("Annotated proteins!")

    if prot_table_path is not None:
        save_table(prot_table_path, df)
    else:
        return df