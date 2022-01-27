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

import pandas as pd
import streamlit as st

from ..constants.conf import sw1_color, sw2_color, sw1_name, sw2_name, loop_resid_dict
from ..constants.pharm import (
    sp2_name,
    none_pharm_name,
    other_pharm_name,
    pharm_color_dict,
)
from ..constants.prot import other_prot_name, none_prot_name, prot_color_dict
from ..functions.lig import lig_lst_dict
from ..functions.table import extract_int, lst_col
from ..functions.lst import str_to_lst
from ..functions.gui import (
    load_st_table,
    show_st_table,
    mask_st_table,
    show_st_structure,
)
from ..functions.col import (
    rename_col_dict,
    pdb_code_col,
    chainid_col,
    bio_lig_col,
    ion_lig_col,
    pharm_lig_col,
    chem_lig_col,
    mod_lig_col,
    mem_lig_col,
    gene_class_col,
    method_col,
    resolution_col,
    r_factor_col,
    space_col,
    mut_status_col,
    nuc_class_col,
    prot_class_col,
    pharm_class_col,
    match_class_col,
    pocket_class_col,
    interf_class_col,
    bound_prot_col,
    bound_prot_chainid_col,
    sw1_col,
    sw2_col,
)


def pdb_page():

    df = load_st_table(__file__)

    st.sidebar.markdown("## PDB Selection")

    pdb_code = st.sidebar.selectbox(
        "Select Entry", [x.upper() for x in lst_col(df, pdb_code_col, unique=True)]
    )
    pdb_df = mask_st_table(df, {pdb_code_col: pdb_code.lower()})

    chainid = st.sidebar.selectbox("Select Chain", lst_col(pdb_df, chainid_col))
    chainid_df = mask_st_table(pdb_df, {chainid_col: chainid})

    st.markdown(
        f"### PDB: [{pdb_code}](https://www.rcsb.org/structure/{pdb_code}) ({df.at[0,gene_class_col]}) - Chain {chainid}"
    )
    st.markdown("---")

    left_col, right_col = st.columns(2)

    left_col.markdown("#### General Information")
    for col in [method_col, resolution_col, r_factor_col, space_col]:
        left_col.markdown(f"**{rename_col_dict[col]}:** {chainid_df.at[0,col]}")

    right_col.markdown("#### Molecular Annotations")

    annot_df = pd.DataFrame()
    for i, col in enumerate(
        [
            mut_status_col,
            nuc_class_col,
            prot_class_col,
            pharm_class_col,
            match_class_col,
            pocket_class_col,
            interf_class_col,
        ]
    ):
        annot_df.at[i, "Molecular Content"] = rename_col_dict[col]
        annot_df.at[i, "Annotation"] = chainid_df.at[0, col]

    show_st_table(annot_df, st_col=right_col)

    for col in [bound_prot_col, bound_prot_chainid_col]:
        right_col.markdown(f"**{rename_col_dict[col]}:** {chainid_df.at[0,col]}")

    sw1_conf = chainid_df.at[0, sw1_col]
    sw2_conf = chainid_df.at[0, sw2_col]

    sw1_html = f'<p style="font-family:sans-serif; color:{sw1_color}; font-size: large;">{sw1_conf}</p>'
    sw2_html = f'<p style="font-family:sans-serif; color:{sw2_color}; font-size: large;">{sw2_conf}</p>'

    left_col.markdown("#### SW1 Conformation")
    left_col.markdown(sw1_html, unsafe_allow_html=True)
    left_col.markdown("#### SW2 Conformation")
    left_col.markdown(sw2_html, unsafe_allow_html=True)

    st.markdown("---")

    left_view_col, right_view_col = st.columns(2)

    style_dict = {"Ribbon": "oval", "Trace": "trace"}

    cartoon_style = style_dict[
        right_view_col.radio("Cartoon Style", ["Ribbon", "Trace"])
    ]

    surf_trans = right_view_col.slider(
        "Surface Transparency", min_value=0.0, max_value=1.0
    )

    left_view_col.markdown("#### Bound Ligand(s)")

    lig_check_dict = dict()
    for col in [
        bio_lig_col,
        ion_lig_col,
        pharm_lig_col,
        chem_lig_col,
        mod_lig_col,
        mem_lig_col,
    ]:
        lig_lst = str_to_lst(chainid_df.at[0, col])
        if "None" not in lig_lst:
            lig_check_dict[col] = dict()
            for lig in lig_lst:
                lig_check_dict[col][lig] = left_view_col.checkbox(
                    f"{rename_col_dict[col]}: {lig}"
                )

    if len(lig_check_dict.keys()) == 0:
        left_view_col.write("No bound ligands.")

    left_view_col.markdown("#### Mutation Site(s)")

    mut_check_dict = dict()
    for mut in str_to_lst(chainid_df.at[0, mut_status_col]):
        if mut != "WT":
            mut_check_dict[mut] = left_view_col.checkbox(mut)

    if len(mut_check_dict.keys()) == 0:
        left_view_col.write("Not Mutated.")

    style_lst = list()
    label_lst = list()

    rotate_view = right_view_col.checkbox("Rotate Structure")

    opacity = 0
    if len(pdb_df) > 1:
        all_chains = right_view_col.checkbox("Show All Chains")

        if all_chains:
            opacity = 0.5

    style_lst.append(
        [
            {
                "chain": chainid,
                "invert": True,
            },
            {
                "cartoon": {
                    "color": "white",
                    "style": cartoon_style,
                    "thickness": 0.2,
                    "opacity": 1,
                    "opacity": opacity,
                }
            },
        ]
    )

    for mut in list(mut_check_dict.keys()):
        style_lst.append(
            [
                {"chain": chainid, "resi": [extract_int(mut)], "atom": "CA"},
                {"sphere": {"color": "red", "radius": 0.8}},
            ]
        )
        if mut in list(mut_check_dict.keys()):
            if mut_check_dict[mut]:
                label_lst.append(
                    [
                        mut,
                        {
                            "backgroundColor": "lightgray",
                            "fontColor": "black",
                            "backgroundOpacity": 0.5,
                        },
                        {"chain": chainid, "resi": [extract_int(mut)], "atom": "CA"},
                    ]
                )

    for lig_col, lig_lst in lig_lst_dict.items():
        lig_lst = str_to_lst(chainid_df.at[0, lig_col])
        if "None" not in lig_lst:
            lig_style = "stick"
            lig_color = "whiteCarbon"
            lig_scheme = "colorscheme"
            lig_radius = "0.2"
            if lig_col == ion_lig_col:
                lig_style = "sphere"
                lig_scheme = "color"
                lig_color = "chartreuse"
                lig_radius = 0.8
            for lig in lig_lst:
                if lig_col != pharm_lig_col:
                    if lig_col == mem_lig_col:
                        style_lst.append(
                            [
                                {
                                    "resn": [lig],
                                },
                                {
                                    lig_style: {
                                        lig_scheme: lig_color,
                                        "radius": lig_radius,
                                    }
                                },
                            ]
                        )

                    else:
                        style_lst.append(
                            [
                                {
                                    "chain": chainid,
                                    "resn": [lig],
                                },
                                {
                                    lig_style: {
                                        lig_scheme: lig_color,
                                        "radius": lig_radius,
                                    }
                                },
                            ]
                        )

                if lig_col in list(lig_check_dict.keys()):
                    if lig_check_dict[lig_col][lig]:
                        label_lst.append(
                            [
                                lig,
                                {
                                    "backgroundColor": "lightgray",
                                    "fontColor": "black",
                                    "backgroundOpacity": 0.5,
                                },
                                {
                                    "chain": chainid,
                                    "resn": lig,
                                },
                            ]
                        )

    pharm_class_lst = str_to_lst(chainid_df.at[0, pharm_class_col])

    if none_pharm_name not in pharm_class_lst:
        if len(pharm_class_lst) > 1:
            pharm_class = other_pharm_name
        else:
            pharm_class = pharm_class_lst[0]

        pharm_color = pharm_color_dict[pharm_class]

        pharm_lig_lst = str_to_lst(chainid_df.at[0, pharm_lig_col])

        for pharm_lig in pharm_lig_lst:
            style_lst.append(
                [
                    {
                        "chain": chainid,
                        "resn": [pharm_lig],
                        "elem": "C",
                    },
                    {"stick": {"color": pharm_color, "radius": 0.2}},
                ]
            )
            style_lst.append(
                [
                    {
                        "chain": chainid,
                        "resn": [pharm_lig],
                        "elem": ["N", "O", "H"],
                    },
                    {"stick": {"colorscheme": "Carbon", "radius": 0.2}},
                ]
            )

            if pharm_class == sp2_name:
                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "resi": 12,
                        },
                        {"stick": {"colorscheme": "lightgrayCarbon", "radius": 0.2}},
                    ]
                )

    prot_class_lst = str_to_lst(chainid_df.at[0, prot_class_col])

    if none_prot_name not in prot_class_lst:
        if len(prot_class_lst) > 1:
            prot_class = other_prot_name
        else:
            prot_class = prot_class_lst[0]
        prot_color = prot_color_dict[prot_class]

        bound_chainid_lst = str_to_lst(chainid_df.at[0, bound_prot_chainid_col])

        for bound_prot_chainid in bound_chainid_lst:
            style_lst.append(
                [
                    {
                        "chain": bound_prot_chainid,
                    },
                    {
                        "cartoon": {
                            "style": cartoon_style,
                            "color": prot_color,
                            "thickness": 0.2,
                            "opacity": 1,
                        }
                    },
                ]
            )

            prot_label = prot_class
            if len(bound_chainid_lst) > 1:
                prot_label += f" - Chain {bound_prot_chainid}"

            label_lst.append(
                [
                    prot_label,
                    {
                        "backgroundColor": "lightgray",
                        "fontColor": "black",
                        "backgroundOpacity": 0.5,
                    },
                    {"chain": bound_prot_chainid},
                ]
            )

    surface_lst = [
        [
            {"opacity": surf_trans, "color": "white"},
            {"chain": chainid, "hetflag": False},
        ]
    ]

    for loop_name, loop_resids in loop_resid_dict.items():

        if loop_name == sw1_name:
            loop_color = sw1_color
            stick_resid = 32
        elif loop_name == sw2_name:
            loop_color = sw2_color
            stick_resid = 71

        surface_lst.append(
            [
                {"opacity": surf_trans, "color": loop_color},
                {"chain": chainid, "resi": loop_resids, "hetflag": False},
            ]
        )

        style_lst.append(
            [
                {
                    "chain": chainid,
                    "resi": [loop_resids],
                },
                {
                    "cartoon": {
                        "style": cartoon_style,
                        "color": loop_color,
                        "thickness": 0.2,
                    }
                },
            ]
        )

        style_lst.append(
            [
                {
                    "chain": chainid,
                    "resi": [stick_resid],
                    "elem": "C",
                },
                {"stick": {"color": loop_color, "radius": 0.2}},
            ]
        )

        style_lst.append(
            [
                {"chain": chainid, "resi": [stick_resid], "elem": ["O", "N", "H"]},
                {"stick": {"colorscheme": "Carbon", "radius": 0.2}},
            ]
        )

    show_st_structure(
        pdb_code,
        style_lst=style_lst,
        surface_lst=surface_lst,
        label_lst=label_lst,
        cartoon_style=cartoon_style,
        spin_on=rotate_view,
        zoom_dict={"chain": chainid},
    )

    st.markdown("---")
    st.write("Copyright (c) 2022 Mitchell Isaac Parker")
