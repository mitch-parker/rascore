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

import math
from collections import defaultdict
from Bio.PDB import calc_dihedral


from .table import extract_int, fix_val
from .lst import lst_nums
from .coord import resid_to_tuple, has_resid, has_atomid, get_atom_vect, calc_norm_dist, get_resname
from .col import (
    phi_col,
    psi_col,
    omega_col,
    chi1_col,
    chi2_col,
    chi3_col,
    chi4_col,
    chi5_col,
    altchi1_col,
    altchi2_col,
)

p_rot = 'p'
m_rot = 'm'
t_rot = 't'
n_rot = '-'

bb_dict = {
    phi_col: [["C", "N", "CA", "C"], [-1, 0, 0, 0]],
    psi_col: [["N", "CA", "C", "N"], [0, 0, 0, 1]],
    omega_col: [["CA", "C", "N", "CA"], [-1, -1, 0, 0]],
}

sc_dict = {
    chi1_col: {
        "ARG": ["N", "CA", "CB", "CG"],
        "ASN": ["N", "CA", "CB", "CG"],
        "ASP": ["N", "CA", "CB", "CG"],
        "CYS": ["N", "CA", "CB", "SG"],
        "GLN": ["N", "CA", "CB", "CG"],
        "GLU": ["N", "CA", "CB", "CG"],
        "HIS": ["N", "CA", "CB", "CG"],
        "ILE": ["N", "CA", "CB", "CG1"],
        "LEU": ["N", "CA", "CB", "CG"],
        "LYS": ["N", "CA", "CB", "CG"],
        "MET": ["N", "CA", "CB", "CG"],
        "PHE": ["N", "CA", "CB", "CG"],
        "PRO": ["N", "CA", "CB", "CG"],
        "SER": ["N", "CA", "CB", "OG"],
        "THR": ["N", "CA", "CB", "OG1"],
        "TRP": ["N", "CA", "CB", "CG"],
        "TYR": ["N", "CA", "CB", "CG"],
        "VAL": ["N", "CA", "CB", "CG1"],
    },
    altchi1_col: {"VAL": ["N", "CA", "CB", "CG2"]},
    chi2_col: {
        "ARG": ["CA", "CB", "CG", "CD"],
        "ASN": ["CA", "CB", "CG", "OD1"],
        "ASP": ["CA", "CB", "CG", "OD1"],
        "GLN": ["CA", "CB", "CG", "CD"],
        "GLU": ["CA", "CB", "CG", "CD"],
        "HIS": ["CA", "CB", "CG", "ND1"],
        "ILE": ["CA", "CB", "CG1", "CD1"],
        "LEU": ["CA", "CB", "CG", "CD1"],
        "LYS": ["CA", "CB", "CG", "CD"],
        "MET": ["CA", "CB", "CG", "SD"],
        "PHE": ["CA", "CB", "CG", "CD1"],
        "PRO": ["CA", "CB", "CG", "CD"],
        "TRP": ["CA", "CB", "CG", "CD1"],
        "TYR": ["CA", "CB", "CG", "CD1"],
    },
    altchi2_col: {
        "ASP": ["CA", "CB", "CG", "OD2"],
        "LEU": ["CA", "CB", "CG", "CD2"],
        "PHE": ["CA", "CB", "CG", "CD2"],
        "TYR": ["CA", "CB", "CG", "CD2"],
    },
    chi3_col: {
        "ARG": ["CB", "CG", "CD", "NE"],
        "GLN": ["CB", "CG", "CD", "OE1"],
        "GLU": ["CB", "CG", "CD", "OE1"],
        "LYS": ["CB", "CG", "CD", "CE"],
        "MET": ["CB", "CG", "SD", "CE"],
    },
    chi4_col: {
        "ARG": ["CG", "CD", "NE", "CZ"],
        "LYS": ["CG", "CD", "CE", "NZ"],
    },
    chi5_col: {"ARG": ["CD", "NE", "CZ", "NH1"]},
}


def get_rama_type(phi, psi, omega):

    phi = float(phi)
    psi = float(psi)
    omega = float(omega)

    omega = math.radians(omega)

    if phi > 180.0:
        phi = -360.0 + phi
    if psi > 180.0:
        psi = -360.0 + psi
    if phi == -180:
        phi = 180
    if psi == -180:
        psi = 180

    regions = defaultdict(list)

    regions["B"].append([-180, -100, 50, 180])
    regions["B"].append([-180, -100, -180, -100])
    regions["A"].append([-180, -100, -100, 50])
    regions["B"].append([-100, 0, -180, -100])
    regions["B"].append([-100, 0, 50, 180])
    regions["A"].append([-100, 0, -100, 50])
    regions["E"].append([0, 150, -180, -50])
    regions["E"].append([0, 150, 100, 180])
    regions["B"].append([150, 180, 100, 180])
    regions["B"].append([150, 180, -180, -50])
    regions["L"].append([0, 180, -50, 100])

    for region in regions:

        for rama in regions[region]:

            if (rama[0] < phi <= rama[1]) and (
                rama[2] < psi <= rama[3] and math.cos(omega) <= 0
            ):

                rama_type = region

            if (rama[0] < phi <= rama[1]) and (
                rama[2] < psi <= rama[3] and math.cos(omega) > 0
            ):

                region = region.lower()
                rama_type = region

    return rama_type


def get_rot_type(rot):

    rot = float(rot)

    if rot == 999.00:
        return n_rot
    elif 0 <= rot < 120:
        return p_rot
    elif 120 <= rot < 240:
        return t_rot
    elif 240 <= rot < 360:
        return m_rot 


def calc_dih_angle(
    structure,
    chainid,
    resid_1,
    resid_2,
    resid_3,
    resid_4,
    atomid_1,
    atomid_2,
    atomid_3,
    atomid_4,
    modelid=None,
    check_dist=True,
    check_resid=True,
    check_atomid=True,
    max_dist=2,
):

    if modelid is None:
        modelid = 0

    index_lst = lst_nums(0, 3)

    resid_lst = [resid_1, resid_2, resid_3, resid_4]
    atomid_lst = [atomid_1, atomid_2, atomid_3, atomid_4]

    for i, resid in enumerate(resid_lst):
        resid_lst[i] = resid_to_tuple(resid)

    structure_angle = 999.00

    calc_vect = True

    for index in index_lst:

        resid = resid_lst[index]
        atomid = atomid_lst[index]

        resid_status = True
        if check_resid:
            resid_status = has_resid(structure, chainid, resid, modelid=modelid)

        atomid_status = True
        if resid_status:
            if check_atomid:
                atomid_status = has_atomid(
                    structure, chainid, resid, atomid, modelid=modelid
                )

        if not resid_status or not atomid_status:
            calc_vect = False
            break

    if calc_vect:

        vect_lst = list()

        for index in index_lst:

            resid = resid_lst[index]
            atomid = atomid_lst[index]

            vect = get_atom_vect(
                structure=structure,
                modelid=modelid,
                chainid=chainid,
                resid=resid,
                atomid=atomid,
            )

            vect_lst.append(vect)

        calc_angle = True

        if check_dist:
            for index in index_lst:
                if index != 0 and index != 3:
                    atom_dist = calc_norm_dist(vect_lst[index], vect_lst[index + 1])
                    if atom_dist > max_dist:
                        calc_angle = False

        if calc_angle:
            structure_angle = round(
                math.degrees(
                    calc_dihedral(vect_lst[0], vect_lst[1], vect_lst[2], vect_lst[3])
                ),
                2,
            )

    return structure_angle


def calc_bb_angle(
    structure,
    chainid,
    curr_resid,
    angle,
    modelid=None,
    prev_resid=None,
    next_resid=None,
    check_dist=True,
    check_resid=True,
    check_atomid=True,
    max_dist=2,
):

    atomid_list = bb_dict[angle][0]
    add_list = bb_dict[angle][1]

    index_lst = lst_nums(0, 3)

    if type(curr_resid) == tuple:
        curr_resnum = curr_resid[1]
    else:
        curr_resnum = extract_int(curr_resid)
        curr_resid = resid_to_tuple(curr_resid)

    resid_lst = list()

    check_dist = False

    for index in index_lst:

        add = add_list[index]

        resid = curr_resid

        if add == -1:
            if prev_resid is not None:
                resid = prev_resid
                if type(prev_resid) == tuple:
                    prev_resnum = prev_resid[1]
                if curr_resnum - prev_resnum != 1:
                    check_dist = True
            else:
                resid = resid_to_tuple(curr_resid + add)
        elif add == 1:
            if next_resid is not None:
                resid = next_resid
                if type(next_resid) == tuple:
                    next_resnum = next_resid[1]
                if next_resnum - curr_resnum != 1:
                    check_dist = True
            else:
                resid = resid_to_tuple(curr_resid + add)

        resid_lst.append(resid)

    bb_angle = calc_dih_angle(
        structure=structure,
        modelid=modelid,
        chainid=chainid,
        resid_1=resid_lst[0],
        resid_2=resid_lst[1],
        resid_3=resid_lst[2],
        resid_4=resid_lst[3],
        atomid_1=atomid_list[0],
        atomid_2=atomid_list[1],
        atomid_3=atomid_list[2],
        atomid_4=atomid_list[3],
        check_dist=check_dist,
        check_resid=check_resid,
        check_atomid=check_atomid,
        max_dist=max_dist,
    )

    return bb_angle


def calc_sc_angle(
    structure,
    chainid,
    resid,
    angle,
    modelid=None,
    check_dist=True,
    check_resid=True,
    check_atomid=True,
    max_dist=2,
):

    if modelid is None:
        modelid = 0

    sc_angle = 999.00

    resid = resid_to_tuple(resid)

    if has_resid(structure, chainid, resid, modelid=modelid):

        residue = structure[fix_val(modelid, return_int=True)][chainid][resid]

        resname = get_resname(residue, letter=False)

        sc_angle_dict = sc_dict[angle]

        if resname in list(sc_angle_dict.keys()):

            atomid_lst = sc_dict[angle][resname]

            sc_angle = calc_dih_angle(
                structure=structure,
                modelid=modelid,
                chainid=chainid,
                resid_1=resid,
                resid_2=resid,
                resid_3=resid,
                resid_4=resid,
                atomid_1=atomid_lst[0],
                atomid_2=atomid_lst[1],
                atomid_3=atomid_lst[2],
                atomid_4=atomid_lst[3],
                check_dist=check_dist,
                check_resid=check_resid,
                check_atomid=check_atomid,
                max_dist=max_dist,
            )

            if sc_angle != 999.00:
                if (angle == chi2_col and resname in ["TYR", "PHE", "ASP"]) or (
                    angle == chi3_col and resname == "GLU"
                ):
                    if sc_angle < 0.00:
                        sc_angle += 180.00
                    elif sc_angle > 180.00:
                        sc_angle -= 180.00
                else:
                    if sc_angle < 0.00:
                        sc_angle += 360.00

    return sc_angle
