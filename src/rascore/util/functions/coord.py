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

import gzip
import numpy as np

from Bio.PDB import (
    MMCIFParser,
    PDBParser,
    MMCIF2Dict,
    MMCIFIO,
    PDBIO,
    Selection,
    NeighborSearch,
    Superimposer,
    HSExposure,
    calc_angle,
)

from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

import math

from .seq import pair_seq_aln
from .table import (
    fix_val,
    replace_str,
    extract_int,
    extract_str,
    make_dict,
    is_int,
    is_str,
)
from .lst import (
    str_to_lst,
    res_to_lst,
    type_lst,
    sort_lst,
    lst_nums,
    lst_unique,
    lst_inter,
)


aa_lst = [
    "CYS",
    "ASP",
    "SER",
    "GLN",
    "LYS",
    "ILE",
    "PRO",
    "THR",
    "PHE",
    "ASN",
    "GLY",
    "HIS",
    "LEU",
    "ARG",
    "TRP",
    "ALA",
    "VAL",
    "GLU",
    "TYR",
    "MET",
]

bb_hb_atomid_dict = dict()
for aa in aa_lst:
    bb_hb_atomid_dict[aa] = ["O", "N"]

h_bb_hb_atomid_dict = dict()
for aa in aa_lst:
    h_bb_hb_atomid_dict[aa] = ["O", "H"]

sc_hb_atomid_dict = {
    "ASP": ["OD1"],
    "GLU": ["OE1"],
    "ASN": ["OD1"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "LYS": ["NZ"],
    "ARG": ["NE", "NH1", "NH2"],
    "GLU": ["OE1", "OE2"],
    "GLN": ["OE1", "NE2"],
    "ASN": ["ND2"],
    "HIS": ["ND1", "NE2"],
}

sc_hb_atomid_lst = list(sc_hb_atomid_dict.keys())


h_sc_hb_atomid_dict = {
    "ASP": ["OD1"],
    "GLU": ["OE1"],
    "ASN": ["OD1", "HD21", "HD22"],
    "SER": ["HG"],
    "THR": ["HG1"],
    "TYR": ["HH"],
    "LYS": ["HZ1", "HZ2", "HZ3"],
    "ARG": ["HE", "HH11", "HH12" "HH21", "HH22"],
    "GLU": ["OE1", "OE2"],
    "GLN": ["OE1", "HE21", "HE22"],
    "HIS": ["ND1", "NE2", "HD1", "HE2"],
}

hb_name = "H-Bond"
wmhb_name = "WMH-Bond"
no_hb_name = "No Bond"

disorder_dict = {2: True, 1: True, 0: False}


def get_pdb_id(pdb_code, chainid):

    return f"{pdb_code}{chainid}"


def get_pdb_code(pdb_id):

    if len(pdb_id) == 5:
        pdb_id = pdb_id[0:4]

    return pdb_id


def get_pdb_chainid(pdb_id):

    if len(pdb_id) == 5:
        pdb_id = pdb_id[4:5]

    return pdb_id


def build_pdb_code_lst(pdb_id_lst):

    pdb_code_lst = list()
    for pdb_id in pdb_id_lst:
        pdb_code = get_pdb_code(pdb_id)
        if pdb_code not in pdb_code_lst:
            pdb_code_lst.append(pdb_code)

    return pdb_code_lst


def build_pdb_dict(pdb_id_lst):

    pdb_dict = {}

    for pdb_id in pdb_id_lst:

        pdb_code = get_pdb_code(pdb_id)
        chainid = get_pdb_chainid(pdb_id)

        if pdb_code not in pdb_dict.keys():
            pdb_dict[pdb_code] = str_to_lst(chainid)
        else:
            pdb_dict[pdb_code].append(chainid)

    return pdb_dict


def load_coord(path):

    if ".cif" in path:
        parser = MMCIFParser(QUIET=True)

    elif ".pdb" in path:
        parser = PDBParser(QUIET=True)

    if ".gz" in path:
        path = gzip.open(path, "rt")

    return parser.get_structure("PARSER", path)


def load_cif_dict(path):

    if ".gz" in path:
        path = gzip.open(path, "rt")

    return MMCIF2Dict.MMCIF2Dict(path)


def search_cif_dict(cif_dict, key):

    return replace_str(str(cif_dict[key][0]), ["?", "-1"], replace=None)


def save_coord(path, structure, sele=None):

    if ".cif" in path:

        io = MMCIFIO()

    elif ".pdb" in path:

        io = PDBIO()

    io.set_structure(structure)

    if sele is not None:
        io.save(path, sele)
    else:
        io.save(path)


def resname_to_letter(resname):

    return aa3to1.get(resname, "X")


def get_resname(residue, letter=False):

    resname = str(residue.get_resname())

    resname = resname.strip()

    if letter:
        resname = resname_to_letter(resname)

    return resname


def unfold_structure(structure, level="A"):

    return Selection.unfold_entities(structure, level)


def get_models(structure):

    return unfold_structure(structure, level="M")


def get_chains(structure):

    return unfold_structure(structure, level="C")


def get_residues(structure):

    return unfold_structure(structure, level="R")


def get_atoms(structure):

    return unfold_structure(structure, level="A")


def get_parent(structure):

    return structure.get_parent()


def get_modelid(model):

    return int(model.id)


def get_chainid(chain):

    return str(chain.id)


def get_resid(residue):

    return residue.id


def get_resmodelid(residue):

    return int(get_modelid(get_parent(get_parent(residue))))


def get_reschainid(residue):

    return str(get_chainid(get_parent(residue)))


def get_resnum(residue):

    return int(get_resid(residue)[1])


def get_resid_str(residue):

    resid = get_resid(residue)

    return f"{resid[1]}{resid[2]}".replace(" ", "")


def resid_to_tuple(resid):

    if type(resid) != tuple:
        resnum = extract_int(resid)
        insert = extract_str(resid)
        if insert == "":
            insert = " "
        resid = tuple([" ", resnum, insert])

    return resid


def get_atomid(atom):

    return str(atom.id)


def is_aa(residue):

    return get_resid(residue)[0] == " "


def is_wat(residue):

    return get_resid(residue)[0] == "W"


def is_het(residue):

    return not is_aa(residue) and not is_wat(residue)


def is_disordered(residue):

    return disorder_dict[residue.is_disordered()]


def has_modelid(structure, modelid):

    if modelid is None:
        modelid = 0

    return structure.has_id(fix_val(modelid, return_int=True))


def has_chainid(structure, chainid, modelid=None):

    if modelid is None:
        modelid = 0

    chainid_status = False
    if has_modelid(structure, modelid):
        chainid_status = structure[fix_val(modelid, return_int=True)].has_id(chainid)

    return chainid_status


def has_resid(structure, chainid, resid, modelid=None):

    if modelid is None:
        modelid = 0

    resid_status = False
    if has_chainid(structure, chainid, modelid=modelid):
        resid_status = structure[fix_val(modelid, return_int=True)][chainid].has_id(
            resid
        )
    return resid_status


def has_atomid(structure, chainid, resid, atomid, modelid=None):

    if modelid is None:
        modelid = 0

    atomid_status = False
    if has_resid(structure, chainid, resid, modelid=modelid):
        atomid_status = structure[fix_val(modelid, return_int=True)][chainid][
            resid
        ].has_id(atomid)
    return atomid_status


def has_altloc(structure, chainid, resid, atomid, altloc, modelid=None):

    if modelid is None:
        modelid = 0

    altloc_status = False

    if has_atomid(structure, chainid, resid, atomid, modelid=modelid):
        atom = structure[fix_val(modelid, return_int=True)][chainid][resid][atomid]
        if is_disordered(atom):
            altloc_status = atom.disordered_has_id(altloc)

    return altloc_status


def get_neighbors(structure):

    return NeighborSearch(get_atoms(structure))


def get_chain_cont(neighbors, chain, max_dist=5, level="A", remove_disordered=False):

    atom_lst = list()

    for residue in chain:

        atom_lst += get_residue_cont(
            neighbors,
            residue,
            max_dist=max_dist,
            level=level,
            remove_disordered=remove_disordered,
        )

    return atom_lst


def get_residue_cont(
    neighbors, residue, max_dist=5, level="A", remove_disordered=False
):

    atom_lst = list()

    for atom in residue:

        atom_lst += get_atom_cont(
            neighbors,
            atom,
            max_dist=max_dist,
            level=level,
            remove_disordered=remove_disordered,
        )

    return atom_lst


def get_atom_cont(neighbors, atom, max_dist=5, level="A", remove_disordered=False):

    atom_lst = neighbors.search(atom.get_coord(), max_dist, level=level)

    if remove_disordered:
        return [x for x in atom_lst if not is_disordered(x)]
    else:
        return atom_lst


def get_seq_lst(structure):

    aa_info = lambda r: (r.id[1], resname_to_letter(get_resname(r)))
    return [aa_info(r) for r in get_residues(structure) if is_aa(r)]


def join_seq_lst(seq_lst):

    return "".join([i[1] for i in seq_lst])


def map_seqs(
    ref_seq_lst,
    mob_seq_lst,
):

    ref_seq = join_seq_lst(ref_seq_lst)
    mob_seq = join_seq_lst(mob_seq_lst)

    aln = pair_seq_aln(ref_seq, mob_seq)

    best_aln = aln[0]

    aln_ref = best_aln[0]
    aln_mob = best_aln[1]

    map_dict = {}
    aa_i_ref, aa_i_mob = 0, 0
    for aa_aln_ref, aa_aln_mob in zip(aln_ref, aln_mob):
        if aa_aln_ref == "-":
            if aa_aln_mob != "-":
                aa_i_mob += 1
        elif aa_aln_mob == "-":
            if aa_aln_ref != "-":
                aa_i_ref += 1
        else:
            assert ref_seq_lst[aa_i_ref][1] == aa_aln_ref
            assert mob_seq_lst[aa_i_mob][1] == aa_aln_mob
            map_dict[ref_seq_lst[aa_i_ref][0]] = mob_seq_lst[aa_i_mob][0]
            aa_i_ref += 1
            aa_i_mob += 1

    return map_dict


def build_map_dict(
    ref_chain,
    mob_chain,
    ref_seq_lst=None,
    mob_seq_lst=None,
    pair_aln=False,
):

    if pair_aln:
        if ref_seq_lst == None:
            ref_seq_lst = get_seq_lst(ref_chain)

        if mob_seq_lst == None:
            mob_seq_lst = get_seq_lst(mob_chain)

        map_dict = map_seqs(
            ref_seq_lst,
            mob_seq_lst,
        )
    else:
        ref_resid_lst = [get_resnum(residue) for residue in ref_chain if is_aa(residue)]
        mob_resid_lst = [get_resnum(residue) for residue in mob_chain if is_aa(residue)]

        resid_lst = [x for x in ref_resid_lst if x in mob_resid_lst]

        map_dict = make_dict(resid_lst, resid_lst)

    return map_dict


def remap_dict(map_dict, resids):

    resid_lst = res_to_lst(resids)

    remap_dict = dict()

    for ref_resid, mob_resid in map_dict.items():
        if ref_resid in resid_lst:
            remap_dict[ref_resid] = mob_resid

    return remap_dict


def get_sup_atoms(ref_chain, mob_chain, map_dict, sup_resids=None):

    if sup_resids is not None:
        map_dict = remap_dict(map_dict, sup_resids)

    ref_ca_lst = list()
    mob_ca_lst = list()

    for ref_resid in map_dict:
        if ref_chain[ref_resid].has_id("CA") and mob_chain[map_dict[ref_resid]].has_id(
            "CA"
        ):
            ref_ca_lst.append(ref_chain[ref_resid]["CA"])
            mob_ca_lst.append(mob_chain[map_dict[ref_resid]]["CA"])

    return (ref_ca_lst, mob_ca_lst)


def sup_coord(mob_structure, ref_ca_lst, mob_ca_lst):

    sup = Superimposer()
    sup.set_atoms(ref_ca_lst, mob_ca_lst)
    sup.apply(mob_structure.get_atoms())

    return mob_structure


def sup_with_map(ref_chain, mob_chain, map_dict, sup_resids=None):

    sup_atoms = get_sup_atoms(ref_chain, mob_chain, map_dict, sup_resids=sup_resids)

    sup_chain = sup_coord(mob_chain, sup_atoms[0], sup_atoms[1])

    return (ref_chain, sup_chain)


def sup_without_map(
    ref_chain,
    mob_chain,
    ref_seq_lst=None,
    mob_seq_lst=None,
    pair_aln=False,
    sup_resids=None,
):

    map_dict = build_map_dict(
        ref_chain,
        mob_chain,
        ref_seq_lst=ref_seq_lst,
        mob_seq_lst=mob_seq_lst,
        pair_aln=pair_aln,
    )

    ref_chain, sup_chain = sup_with_map(
        ref_chain, mob_chain, map_dict, sup_resids=sup_resids
    )

    return (ref_chain, sup_chain, map_dict)


def calc_norm_dist(vect_1, vect_2, dist_dec=2):

    return round(np.linalg.norm(vect_1 - vect_2), dist_dec)


def prep_resid(structure, chainid, resid, modelid=None):

    if modelid is None:
        modelid = 0

    if is_str(resid):
        for residue in structure[fix_val(modelid, return_int=True)][chainid]:
            if is_het(residue):
                if resid == get_resname(residue):
                    resid = get_resid(residue)
    elif is_int(resid):
        resid = resid_to_tuple(resid)

    return resid


def select_atom_altloc(atom, altloc):

    if altloc is not None:
        atom.disordered_select(altloc)

    return atom


def get_atom_vect(
    structure,
    chainid,
    resid,
    atomid,
    modelid=None,
    altloc=None,
    coord=False,
):

    if modelid is None:
        modelid = 0

    atom_vect = 999.00

    resid = prep_resid(structure, chainid, resid, modelid=modelid)

    if has_atomid(structure, chainid, resid, atomid, modelid=modelid):
        atom = structure[fix_val(modelid, return_int=True)][chainid][resid][atomid]
        if altloc is not None:
            if has_altloc(structure, chainid, resid, atomid, altloc, modelid=modelid):
                atom = select_atom_altloc(atom, altloc)
        if coord:
            atom_vect = atom.get_coord()
        else:
            atom_vect = atom.get_vector()

    return atom_vect


def calc_rmsd(
    ref_chain,
    sup_chain,
    map_dict=None,
    rmsd_resids=None,
    rmsd_atomids="CA",
    pair_aln=False,
):

    if map_dict is None:
        map_dict = build_map_dict(ref_chain, sup_chain, pair_aln=pair_aln)

    map_resid_lst = list(map_dict.keys())

    if rmsd_resids is None:
        resid_lst = map_resid_lst
    else:
        resid_lst = res_to_lst(rmsd_resids)

    atomid_lst = type_lst(rmsd_atomids)

    ref_vect_lst = list()
    sup_vect_lst = list()

    for resid in resid_lst:
        if resid in map_resid_lst:
            ref_resid = resid_to_tuple(resid)
            sup_resid = resid_to_tuple(map_dict[resid])
            for atomid in atomid_lst:
                ref_vect = 999.00
                if ref_chain[ref_resid].has_id(atomid):
                    ref_vect = ref_chain[ref_resid][atomid].get_coord()
                if type(ref_vect) != float:
                    sup_vect = 999.00
                    if sup_chain[sup_resid].has_id(atomid):
                        sup_vect = sup_chain[sup_resid][atomid].get_coord()
                    if type(sup_vect) != float:
                        ref_vect_lst.append(list(ref_vect))
                        sup_vect_lst.append(list(sup_vect))

    tot_vects = len(ref_vect_lst)

    if tot_vects == 0:
        rmsd_val = 999.00
    else:
        ref_vects = np.array(ref_vect_lst)
        sup_vects = np.array(sup_vect_lst)
        diff_vect = ref_vects - sup_vects
        rmsd_val = float(np.sqrt((diff_vect ** 2).sum() / tot_vects))

    return rmsd_val


def build_add_resid_lst(
    coord_path,
    modelid,
    chainid,
    resid_lst,
    dict_resid_lst,
    max_ca_dist=4,
    ext_mult=1,
):

    min_resid = resid_lst[0]
    max_resid = resid_lst[-1]

    search_resid_lst = resid_lst.copy()

    search_ext = len(resid_lst * ext_mult)

    search_lst = lst_nums(1, search_ext)

    for i in search_lst:
        search_resid_lst.append(min_resid - i)
        search_resid_lst.append(max_resid + i)

    search_resid_lst = sort_lst(search_resid_lst)

    search_resid_dict = dict()

    for search_resid in search_resid_lst:
        search_resid_dict[search_resid] = [
            x for x in dict_resid_lst if extract_int(x) == extract_int(search_resid)
        ]

    loaded = False
    add_resid_lst = list()
    for resid in resid_lst:
        curr_resid_lst = search_resid_dict[resid]
        if len(curr_resid_lst) == 0:
            prev_resid = None
            next_resid = None
            for i in search_lst:
                prev_resid_lst = search_resid_dict[resid - i]
                next_resid_lst = search_resid_dict[resid + i]
                if prev_resid is None:
                    if len(prev_resid_lst) > 0:
                        prev_resid = prev_resid_lst[-1]
                        curr_resid = extract_int(prev_resid) + i
                if next_resid is None:
                    if len(next_resid_lst) > 0:
                        next_resid = next_resid_lst[0]
                        curr_resid = extract_int(next_resid) - i
                if prev_resid is not None and next_resid is not None:
                    if not loaded:
                        structure = load_coord(coord_path)
                        loaded = True

                    atom_dist = calc_atom_dist(
                        structure=structure,
                        chainid_1=chainid,
                        resid_1=prev_resid,
                        chainid_2=chainid,
                        resid_2=next_resid,
                        atomid_1="CA",
                        atomid_2="CA",
                        modelid_1=modelid,
                        modelid_2=modelid,
                    )
                    if atom_dist > max_ca_dist:
                        curr_resid_lst = type_lst(str(curr_resid))
                    break

        for curr_resid in curr_resid_lst:
            if curr_resid not in add_resid_lst:
                add_resid_lst.append(curr_resid)

    return add_resid_lst


def get_hb_atomid_lst(
    structure, chainid, resid, modelid=None, hb_sc=True, hb_bb=True, use_h=False
):

    if modelid is None:
        modelid = 0

    atomid_lst = list()

    used_bb_hb_atomid_dict = bb_hb_atomid_dict.copy()
    used_sc_hb_atomid_dict = sc_hb_atomid_dict.copy()
    if use_h:
        used_bb_hb_atomid_dict = h_bb_hb_atomid_dict.copy()
        used_sc_hb_atomid_dict = h_sc_hb_atomid_dict.copy()

    if has_resid(structure, chainid, resid, modelid=modelid):

        residue = structure[fix_val(modelid, return_int=True)][chainid][resid]

        resname = get_resname(residue)

        if is_aa(residue):

            atomid_lst = list()

            if hb_bb is not False:
                if resname in aa_lst:
                    atomid_lst += used_bb_hb_atomid_dict[resname]
                    if type(hb_bb) == str:
                        atomid_lst = [x for x in atomid_lst if hb_bb in x]
            if hb_sc is not False:
                if resname in sc_hb_atomid_lst:
                    atomid_lst += used_sc_hb_atomid_dict[resname]
                    if type(hb_sc) == str:
                        atomid_lst = [x for x in atomid_lst if hb_sc in x]
        else:
            for atom in residue:
                atomid = get_atomid(atom)
                if atomid[0] == "O" or atomid[0] == "N" or atomid[0] == "H":
                    atomid_lst.append(atomid)

    if len(atomid_lst) == 0:
        atomid_lst = None

    return atomid_lst


def pred_wmhb_dist(max_wmhb_dist, max_wmhb_angle):

    return np.sqrt(
        (max_wmhb_dist ** 2)
        + (max_wmhb_dist ** 2)
        - (2 * max_wmhb_dist * max_wmhb_dist * math.cos(math.radians(max_wmhb_angle)))
    )


def get_adj_atom(
    structure, chainid, resid, atomid, altloc=None, modelid=None, adj_dist=4
):

    if modelid is None:
        modelid == 0

    adj_atom = None

    if has_resid(structure, chainid, resid, modelid=modelid):
        chain = structure[fix_val(modelid, return_int=True)][chainid]
        atom = chain[resid][atomid]

        neighbors = get_neighbors(chain)
        atom = select_atom_altloc(atom, altloc)

        atom_vect = atom.get_vector()

        adj_dist_lst = list()
        adj_atom_lst = list()

        for adj_atom in get_atom_cont(neighbors, atom, max_dist=adj_dist):
            adj_atomid = get_atomid(adj_atom)
            if adj_atomid != atomid:
                if (atomid[0] == "H") or (
                    atomid[0] != "H"
                    and (adj_atomid[0] == "C" or is_het(atom.get_parent()))
                ):
                    adj_dist_lst.append(
                        calc_norm_dist(atom_vect, adj_atom.get_vector())
                    )
                    adj_atom_lst.append(adj_atom)

        if len(adj_atom_lst) > 0:
            adj_atom = adj_atom_lst[adj_dist_lst.index(min(adj_dist_lst))]

    return adj_atom


def find_hb_status(
    atom_dist,
    structure,
    chainid_1,
    resid_1,
    atomid_1,
    chainid_2,
    resid_2,
    atomid_2,
    altloc_1=None,
    altloc_2=None,
    modelid_1=None,
    modelid_2=None,
    adj_atomid_1=None,
    adj_atomid_2=None,
    min_hb_dist=2.0,
    max_hb_dist=3.2,
    min_wmhb_dist=2.0,
    max_wmhb_dist=3.0,
    min_hb_angle=90,
    max_hb_angle=180,
    min_wmhb_angle=80,
    max_wmhb_angle=140,
    return_angle=False,
):

    if modelid_1 is None:
        modelid_1 = 0
    if modelid_2 is None:
        modelid_2 = 0

    hb_status = no_hb_name
    hb_angle_1 = 999.00
    hb_angle_2 = 999.00
    wmhb_angle = 999.00
    outlier_status = False

    min_wmhb_diag_dist = pred_wmhb_dist(min_wmhb_dist, min_wmhb_angle)
    max_wmhb_diag_dist = pred_wmhb_dist(max_wmhb_dist, max_wmhb_angle)

    if has_resid(structure, chainid_1, resid_1, modelid=modelid_1) and has_resid(
        structure, chainid_2, resid_2, modelid=modelid_2
    ):
        residue_1 = structure[fix_val(modelid_1, return_int=True)][chainid_1][resid_1]
        residue_2 = structure[fix_val(modelid_2, return_int=True)][chainid_2][resid_2]

        resname_1 = get_resname(residue_1)
        resname_2 = get_resname(residue_2)

        calc_hb = True
        if is_aa(residue_1):
            if (atomid_1 not in bb_hb_atomid_dict[resname_1]) and (
                atomid_1 not in h_bb_hb_atomid_dict[resname_1]
            ):
                calc_hb = False
                if resname_1 in sc_hb_atomid_lst:
                    if (atomid_1 in sc_hb_atomid_dict[resname_1]) or (
                        (atomid_1 in h_sc_hb_atomid_dict[resname_1])
                    ):
                        calc_hb = True

        if is_aa(residue_2):
            if (atomid_2 not in bb_hb_atomid_dict[resname_2]) and (
                atomid_2 not in h_bb_hb_atomid_dict[resname_2]
            ):
                calc_hb = False
                if resname_2 in sc_hb_atomid_lst:
                    if (atomid_2 in sc_hb_atomid_dict[resname_2]) or (
                        (atomid_2 in h_sc_hb_atomid_dict[resname_2])
                    ):
                        calc_hb = True

        if calc_hb:

            find_adj_1 = True
            if adj_atomid_1 is not None:
                if has_atomid(
                    structure, chainid_1, resid_1, adj_atomid_1, modelid=modelid_1
                ):
                    adj_atom_1 = structure[fix_val(modelid_1, return_int=True)][
                        chainid_1
                    ][resid_1][adj_atomid_1]
                    find_adj_1 = False
            if find_adj_1:
                adj_atom_1 = get_adj_atom(
                    structure,
                    chainid_1,
                    resid_1,
                    atomid_1,
                    altloc=altloc_1,
                    modelid=modelid_1,
                )

            find_adj_2 = True
            if adj_atomid_2 is not None:
                if has_atomid(
                    structure, chainid_2, resid_2, adj_atomid_2, modelid=modelid_2
                ):
                    adj_atom_2 = structure[fix_val(modelid_2, return_int=True)][
                        chainid_2
                    ][resid_2][adj_atomid_2]
                    find_adj_2 = False
            if find_adj_2:
                adj_atom_2 = get_adj_atom(
                    structure,
                    chainid_2,
                    resid_2,
                    atomid_2,
                    altloc=altloc_2,
                    modelid=modelid_2,
                )

            if adj_atom_1 is not None and adj_atom_2 is not None:

                atom_1 = select_atom_altloc(residue_1[atomid_1], altloc_1)
                atom_2 = select_atom_altloc(residue_2[atomid_2], altloc_2)

                atom_1_vect = atom_1.get_vector()
                atom_2_vect = atom_2.get_vector()
                adj_atom_1_vect = adj_atom_1.get_vector()
                adj_atom_2_vect = adj_atom_2.get_vector()

                if (min_hb_dist <= atom_dist <= max_hb_dist) and (
                    (atomid_1[0] != atomid_2[0])
                    or (resname_1 not in aa_lst or resname_2 not in aa_lst)
                ):

                    hb_angle_1 = math.degrees(
                        calc_angle(adj_atom_1_vect, atom_1_vect, atom_2_vect)
                    )

                    hb_angle_2 = math.degrees(
                        calc_angle(adj_atom_2_vect, atom_2_vect, atom_1_vect)
                    )
                    if (min_hb_angle <= hb_angle_1 <= max_hb_angle) and (
                        min_hb_angle <= hb_angle_2 <= max_hb_angle
                    ):
                        hb_status = hb_name
                    else:
                        outlier_status = True
                else:
                    neighbors_1 = get_neighbors(
                        structure[fix_val(modelid_1, return_int=True)][chainid_1]
                    )
                    neighbors_2 = get_neighbors(
                        structure[fix_val(modelid_2, return_int=True)][chainid_2]
                    )

                    atom_cont_1_lst = get_atom_cont(
                        neighbors_1, atom_1, max_dist=max_wmhb_dist
                    )
                    atom_cont_2_lst = get_atom_cont(
                        neighbors_2, atom_2, max_dist=max_wmhb_dist
                    )

                    residue_cont_1_lst = [x.get_parent() for x in atom_cont_1_lst]

                    residue_cont_2_lst = [x.get_parent() for x in atom_cont_2_lst]

                    wat_cont_lst_1 = lst_unique(
                        [x for x in residue_cont_1_lst if is_wat(x)]
                    )
                    wat_cont_lst_2 = lst_unique(
                        [x for x in residue_cont_2_lst if is_wat(x)]
                    )

                    wat_cont_lst = lst_inter(wat_cont_lst_1, wat_cont_lst_2)

                    if len(wat_cont_lst) > 0:

                        for wat_residue in wat_cont_lst:

                            wat_atom = wat_residue["O"]

                            wat_atom_vect = wat_atom.get_vector()

                            wmhb_dist_1 = calc_norm_dist(atom_1_vect, wat_atom_vect)

                            wmhb_dist_2 = calc_norm_dist(atom_2_vect, wat_atom_vect)

                            if (min_wmhb_dist < wmhb_dist_1 <= max_wmhb_dist) and (
                                min_wmhb_dist < wmhb_dist_2 <= max_wmhb_dist
                            ):

                                wmhb_angle = math.degrees(
                                    calc_angle(
                                        atom_1_vect,
                                        wat_atom_vect,
                                        atom_2_vect,
                                    )
                                )
                                if min_wmhb_angle < wmhb_angle < max_wmhb_angle:
                                    hb_angle_1 = math.degrees(
                                        calc_angle(
                                            adj_atom_1_vect,
                                            atom_1_vect,
                                            wat_atom_vect,
                                        )
                                    )
                                    hb_angle_2 = math.degrees(
                                        calc_angle(
                                            adj_atom_2_vect,
                                            atom_2_vect,
                                            wat_atom_vect,
                                        )
                                    )
                                    if (
                                        min_hb_angle <= hb_angle_1 <= max_hb_angle
                                    ) and (min_hb_angle <= hb_angle_2 <= max_hb_angle):
                                        hb_status = wmhb_name
                                        break
                                    else:
                                        outlier_status = True
                                else:
                                    outlier_status = True

                    elif min_wmhb_diag_dist < atom_dist < max_wmhb_diag_dist:
                        hb_status = wmhb_name

    if return_angle:
        return (hb_status, hb_angle_1, hb_angle_2, wmhb_angle, outlier_status)
    else:
        return hb_status


def get_altloc_lst(structure, chainid, resid, atomid, modelid=None):

    altloc_lst = [None]
    for altloc in ["A", "B", "C", "D", "E"]:
        if has_altloc(structure, chainid, resid, atomid, altloc, modelid=modelid):
            altloc_lst.append(altloc)

    return altloc_lst


def calc_atom_dist(
    structure,
    chainid_1,
    resid_1,
    chainid_2,
    resid_2,
    atomid_1=None,
    atomid_2=None,
    modelid_1=None,
    modelid_2=None,
    adj_atomid_1=None,
    adj_atomid_2=None,
    check_hb=False,
    use_h=False,
    hb_sc=True,
    hb_bb=True,
    min_hb_dist=2.0,
    max_hb_dist=3.2,
    min_wmhb_dist=2.0,
    max_wmhb_dist=3.0,
    min_hb_angle=90,
    max_hb_angle=180,
    min_wmhb_angle=80,
    max_wmhb_angle=140,
    dist_dec=2,
    return_vect=False,
):

    if modelid_1 is None:
        modelid_1 = 0
    if modelid_2 is None:
        modelid_2 = 0

    resid_1 = prep_resid(structure, chainid_1, resid_1, modelid=modelid_1)
    resid_2 = prep_resid(structure, chainid_2, resid_2, modelid=modelid_2)

    if atomid_1 is None:
        atomid_1_lst = get_hb_atomid_lst(
            structure,
            chainid_1,
            resid_1,
            modelid=modelid_1,
            hb_sc=hb_sc,
            hb_bb=hb_bb,
            use_h=use_h,
        )
    else:
        atomid_1_lst = type_lst(atomid_1)

    if atomid_2 is None:
        atomid_2_lst = get_hb_atomid_lst(
            structure,
            chainid_2,
            resid_2,
            modelid=modelid_2,
            hb_sc=hb_sc,
            hb_bb=hb_bb,
            use_h=use_h,
        )
    else:
        atomid_2_lst = type_lst(atomid_2)

    atom_dist_lst = list()
    vect_1_lst = list()
    vect_2_lst = list()

    if atomid_1_lst is not None and atomid_2_lst is not None:

        altloc_1_lst = [
            get_altloc_lst(structure, chainid_1, resid_1, x, modelid=modelid_1)
            for x in atomid_1_lst
        ]
        altloc_2_lst = [
            get_altloc_lst(structure, chainid_2, resid_2, x, modelid=modelid_2)
            for x in atomid_2_lst
        ]

        atomid_pair_lst = list()
        for i, atomid_1 in enumerate(atomid_1_lst):
            for altloc_1 in altloc_1_lst[i]:
                for j, atomid_2 in enumerate(atomid_2_lst):
                    for altloc_2 in altloc_2_lst[j]:
                        atomid_pair_lst.append((atomid_1, altloc_1, atomid_2, altloc_2))

        for atomid_pair in atomid_pair_lst:

            atomid_1 = atomid_pair[0]
            altloc_1 = atomid_pair[1]

            atomid_2 = atomid_pair[2]
            altloc_2 = atomid_pair[3]

            vect_1 = 999.00
            vect_2 = 999.00
            atom_dist = 999.00

            vect_1 = get_atom_vect(
                structure=structure,
                modelid=modelid_1,
                chainid=chainid_1,
                resid=resid_1,
                atomid=atomid_1,
                altloc=altloc_1,
            )

            if vect_1 != 999.00:

                vect_2 = get_atom_vect(
                    structure=structure,
                    modelid=modelid_2,
                    chainid=chainid_2,
                    resid=resid_2,
                    atomid=atomid_2,
                    altloc=altloc_2,
                )

                if vect_2 != 999.00:

                    atom_dist = calc_norm_dist(vect_1, vect_2, dist_dec=dist_dec)

            vect_1_lst.append(vect_1)
            vect_2_lst.append(vect_2)
            atom_dist_lst.append(atom_dist)

    atomid_pair = None
    if len(atom_dist_lst) == 0:
        vect_1 = 999.00
        vect_2 = 999.00
        atom_dist = 999.00
    else:
        atom_dist = min(atom_dist_lst)
        vect_1 = vect_1_lst[atom_dist_lst.index(atom_dist)]
        vect_2 = vect_2_lst[atom_dist_lst.index(atom_dist)]
        if atom_dist != 999.00:
            atomid_pair = atomid_pair_lst[atom_dist_lst.index(atom_dist)]

    if check_hb:

        hb_status = no_hb_name
        hb_angle_1 = 999.00
        hb_angle_2 = 999.00
        wmhb_angle = 999.00
        outlier_status = True

        if atomid_pair is not None:

            atomid_1 = atomid_pair[0]
            altloc_1 = atomid_pair[1]

            atomid_2 = atomid_pair[2]
            altloc_2 = atomid_pair[3]

            hb_result = find_hb_status(
                atom_dist,
                structure,
                chainid_1,
                resid_1,
                atomid_1,
                chainid_2,
                resid_2,
                atomid_2,
                altloc_1=altloc_1,
                altloc_2=altloc_2,
                modelid_1=modelid_1,
                modelid_2=modelid_2,
                adj_atomid_1=adj_atomid_1,
                adj_atomid_2=adj_atomid_2,
                min_hb_dist=min_hb_dist,
                max_hb_dist=max_hb_dist,
                min_wmhb_dist=min_wmhb_dist,
                max_wmhb_dist=max_wmhb_dist,
                min_hb_angle=min_hb_angle,
                max_hb_angle=max_hb_angle,
                min_wmhb_angle=min_wmhb_angle,
                max_wmhb_angle=max_wmhb_angle,
                return_angle=True,
            )

            hb_status = hb_result[0]
            hb_angle_1 = hb_result[1]
            hb_angle_2 = hb_result[2]
            wmhb_angle = hb_result[3]
            outlier_status = hb_result[4]

    result = (atom_dist,)

    if return_vect:
        result += (vect_1, vect_2)

    if check_hb:
        result += (hb_status, hb_angle_1, hb_angle_2, wmhb_angle, outlier_status)

    if len(result) == 1:
        result = result[0]

    return result


def build_hse_dict(model, radius=12):

    return HSExposure.HSExposureCA(model, radius=radius)


def get_hse_val(hse_dict, chainid, resid, index=None):

    hse_tuple = tuple([chainid, resid_to_tuple(resid)])

    hse = 999.00
    if hse_tuple in list(hse_dict.keys()):
        hse = hse_dict[hse_tuple]

    if index is not None:
        hse = hse[index]

    return hse


def get_hse_up(hse_dict, chainid, resid):

    return get_hse_val(hse_dict, chainid, resid, index=0)


def get_hse_down(hse_dict, chainid, resid):

    return get_hse_val(hse_dict, chainid, resid, index=1)
