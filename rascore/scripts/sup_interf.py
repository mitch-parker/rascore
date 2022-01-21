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
from tqdm import tqdm

from ..functions import *


def run_sup(mob_structure, ref_chain, mob_chain, sup_resids=None, pair_aln=False):

    map_dict = build_map_dict(
        ref_chain,
        mob_chain,
        pair_aln=pair_aln,
    )
    sup_atoms = get_sup_atoms(ref_chain, mob_chain, map_dict, sup_resids=sup_resids)

    return sup_coord(mob_structure, sup_atoms[0], sup_atoms[1])


def sup_interf(
    mob_path_lst,
    output_path=None,
    ref_path=None,
    sup_resids=None,
    rmsd_resids=None,
    pair_aln=False,
):

    if output_path is None:
        output_path = os.getcwd()
    output_path += f"/{sup_str}"
    append_path(output_path)

    if ref_path is None:
        ref_path = mob_path_lst[0]

    ref_structure = load_coord(ref_path)

    ref_chain_lst = get_chains(ref_structure)

    sup_path_lst = list()

    for mob_path in tqdm(
        mob_path_lst, desc="Superposing interfaces", position=0, leave=True
    ):
        mob_structure = load_coord(mob_path)
        sup_lst = list()
        rmsd_lst = list()
        for ref_chain in ref_chain_lst:
            ref_chainid = get_chainid(ref_chain)
            for mob_chain in get_chains(mob_structure):
                mob_chainid = get_chainid(mob_chain)
                sup_lst.append((ref_chainid, mob_chainid))
                sup_structure = run_sup(
                    mob_structure,
                    ref_chain,
                    mob_chain,
                    sup_resids=sup_resids,
                    pair_aln=pair_aln,
                )
                chain_rmsd_lst = list()
                for sup_chain in get_chains(sup_structure):
                    chain_rmsd_lst.append(
                        min(
                            [
                                calc_rmsd(
                                    x,
                                    sup_chain,
                                    rmsd_resids=rmsd_resids,
                                )
                                for x in ref_chain_lst
                            ]
                        )
                    )
                rmsd_lst.append(max(chain_rmsd_lst))

        best_sup = sup_lst[rmsd_lst.index(min(rmsd_lst))]

        best_sup_structure = run_sup(
            mob_structure,
            ref_structure[0][best_sup[0]],
            mob_structure[0][best_sup[1]],
            sup_resids=sup_resids,
            pair_aln=pair_aln,
        )

        sup_path = f"{output_path}/{get_file_name(mob_path)}"

        sup_path_lst.append(sup_path)

        save_coord(sup_path, best_sup_structure)

    print("Superposed Interfaces!")

    return sup_path_lst
