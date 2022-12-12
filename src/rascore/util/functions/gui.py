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

import pandas as pd
import uuid
import re
from random import randint
import streamlit as st
import py3Dmol
from stmol import showmol
from io import BytesIO

from .file import entry_table_file
from .table import mask_equal, merge_dicts
from .path import (
    load_table,
    load_json,
    get_file_path,
    get_dir_path,
    get_neighbor_path,
    pages_str,
    data_str,
    functions_str,
)
from .lig import lig_col_lst 
from .lst import res_to_lst, str_to_lst,  type_lst
from .color import red_hex
from .col import (rename_col_dict, date_col, nuc_class_col, match_class_col, prot_class_col, 
                gene_class_col, interf_class_col, pocket_class_col, pdb_code_col, chainid_col, ion_lig_col, 
                bound_prot_chainid_col, pharm_lig_col, mem_lig_col, pharm_class_col)

from ..constants.nuc import nuc_class_lst
from ..constants.pharm import match_class_lst, pocket_class_lst, pharm_color_dict, none_pharm_name, other_pharm_name, sp2_name
from ..constants.prot import prot_class_lst, prot_color_dict, none_prot_name, other_prot_name
from ..constants.conf import sw1_name_lst, sw2_name_lst, y32_name_lst, y71_name_lst, y32_name, y71_name, sw1_name, sw2_name, loop_resid_dict, loop_color_dict
from ..constants.gene import gene_class_lst
from ..constants.dimer import interf_class_lst

class_order_dict = {nuc_class_col: nuc_class_lst, match_class_col: match_class_lst, pocket_class_col: pocket_class_lst,
                    prot_class_col: prot_class_lst, sw1_name: sw1_name_lst, sw2_name: sw2_name_lst,
                    y32_name: y32_name_lst, y71_name: y71_name_lst,
                    gene_class_col: gene_class_lst, interf_class_col: interf_class_lst}

sw1_resid_lst = res_to_lst(loop_resid_dict[sw1_name])
sw2_resid_lst = res_to_lst(loop_resid_dict[sw2_name])

mitch_twitter = '<a href="https://twitter.com/Mitch_P?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">Follow @Mitch_P</a><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>'
roland_twitter = '<a href="https://twitter.com/RolandDunbrack?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">Follow @RolandDunbrack</a><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>'


ribbon_name = "Ribbon"
trace_name = "Trace"
standard_name = "Standard"
aa_name = "Amino Acid"

def write_st_end():

    df = load_table(
        get_file_path(
            entry_table_file,
            dir_path=get_neighbor_path(__file__, functions_str, data_str),
        )
    )

    df[date_col] = pd.to_datetime(df[date_col])
    df[date_col] = df[date_col].dt.strftime("%Y-%m")

    st.markdown("---")
    st.markdown(
        "Developed and Maintained by Mitchell Parker, Bulat Faezov, and Roland Dunbrack"
    )
    st.markdown(
        "[Dunbrack Lab](https://dunbrack.fccc.edu/retro/) - [Fox Chase Cancer Center](https://www.foxchase.org)"
    )
    st.markdown(f"Most Recently Deposited Entry {df[date_col].max()}")
    st.markdown("Copyright (c) 2022 Mitchell Isaac Parker")

def reorder_st_cols(df, row, col):

    class_lst = list(class_order_dict.keys())

    if row in class_lst:
        row_order = [x for x in class_order_dict[row] + ["All"] if x in list(df.index.values)] 

    if col in class_lst:
        col_order = [x for x in class_order_dict[col] + ["All"] if x in list(df.columns)] 

    if row in class_lst and col in class_lst:
        df = df.loc[row_order, col_order]
    elif row not in class_lst and col in class_lst:
        df = df.loc[:, col_order]
    elif row in class_lst and col not in class_lst:
        df = df.loc[row_order, :]
    
    return df

def get_st_file_path(st_file):

    return get_file_path(
        f"{randint(0,3261994)}_{st_file.name}",
        dir_path=get_neighbor_path(__file__, functions_str, data_str),
    )


def save_st_file(st_file):

    st_file_path = get_st_file_path(st_file)
    with open(st_file_path, "wb") as file:
        file.write(st_file.getbuffer())
    return st_file_path

def show_st_fig(fig, st_col=None):

    byt = BytesIO()
    fig.savefig(byt, format="png")
    if st_col is None:
        st.image(byt)
    else:
        st_col.image(byt)


def get_html_text(text_color_dict, font_size="medium", font_weight="normal"):

    html_str = ""
    for text, color in text_color_dict.items():

        size = font_size
        if type(font_size) == dict:
            size = font_size[text]

        weight = font_weight
        if type(font_weight) == dict:
            weight = font_weight[text]

        html_str += f'<span style="font-family:sans-serif; font-size: {size}; font-weight: {weight}; color:{color};">{text}</span>'

    return html_str


def load_st_table(file_path, file_name=None, json_format=False):

    if file_name is None:
        file_name = entry_table_file

    file_path = get_file_path(
            file_name,
            dir_path=get_dir_path(
                dir_path=get_neighbor_path(file_path, pages_str, data_str)
            ))

    if json_format:
        return load_json(file_path)
    else:
        return load_table(file_path)
   
def mask_st_table(df, col_dict):

    mask_df = df.copy()

    for col in list(col_dict.keys()):
        val = type_lst(col_dict[col])
        if len(val) > 0:
            if "All" not in val:
                mask_df = mask_equal(mask_df, col, val)

    return mask_df


def rename_st_cols(df, col_lst=None):

    if col_lst is None:
        col_lst = [x for x in list(rename_col_dict.keys()) if x in list(df.columns)]

    return df.loc[:, col_lst].rename(columns=rename_col_dict)


def show_st_dataframe(df, st_col=None, hide_index=True):

    if hide_index:
        hide_dataframe_row_index = """
                <style>
                .row_heading.level0 {display:none}
                .blank {display:none}
                </style>
                """

        st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)

    if st_col is None:
        st.dataframe(df)
    else:
        st_col.dataframe(df)


def show_st_table(df, st_col=None, hide_index=True):

    if hide_index:
        hide_table_row_index = """
                <style>
                tbody th {display:none}
                .blank {display:none}
                </style>
                """
        st.markdown(hide_table_row_index, unsafe_allow_html=True)

    if st_col is None:
        st.table(df)
    else:
        st_col.table(df)


def create_st_button(link_text, link_url, hover_color="#e78ac3", st_col=None):

    button_uuid = str(uuid.uuid4()).replace("-", "")
    button_id = re.sub("\d+", "", button_uuid)

    button_css = f"""
        <style>
            #{button_id} {{
                background-color: rgb(255, 255, 255);
                color: rgb(38, 39, 48);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial;

            }}
            #{button_id}:hover {{
                border-color: {hover_color};
                color: {hover_color};
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: {hover_color};
                color: white;
                }}
        </style> """

    html_str = f'<a href="{link_url}" target="_blank" id="{button_id}";>{link_text}</a><br></br>'

    if st_col is None:
        st.markdown(button_css + html_str, unsafe_allow_html=True)
    else:
        st_col.markdown(button_css + html_str, unsafe_allow_html=True)


def download_st_file(file_path, file_name, download_text, st_col=None):

    with open(file_path, "rb") as file:
        if st_col is None:
            st.download_button(download_text, file, file_name=file_name)
        else:
            st_col.download_button(download_text, file, file_name=file_name)


def encode_st_df(df):

    return df.to_csv(sep="\t", index=False).encode("utf-8")


def download_st_df(df, file_name, download_text, st_col=None):

    if st_col is None:
        st.download_button(
            label=download_text,
            data=encode_st_df(df),
            file_name=file_name,
        )
    else:
        st_col.download_button(
            label=download_text,
            data=encode_st_df(df),
            file_name=file_name,
        )


def show_st_3dmol(
    pdb_code,
    style_lst=None,
    label_lst=None,
    reslabel_lst=None,
    zoom_dict=None,
    surface_lst=None,
    cartoon_style="trace",
    cartoon_radius=0.2,
    cartoon_color="lightgray",
    zoom=1,
    spin_on=False,
    width=900,
    height=600,
):

    view = py3Dmol.view(query=f"pdb:{pdb_code.lower()}", width=width, height=height)

    view.setStyle(
        {
            "cartoon": {
                "style": cartoon_style,
                "color": cartoon_color,
                "thickness": cartoon_radius,
            }
        }
    )

    if surface_lst is not None:
        for surface in surface_lst:
            view.addSurface(py3Dmol.VDW, surface[0], surface[1])

    if style_lst is not None:
        for style in style_lst:
            view.addStyle(
                style[0],
                style[1],
            )

    if label_lst is not None:
        for label in label_lst:
            view.addLabel(label[0], label[1], label[2])

    if reslabel_lst is not None:
        for reslabel in reslabel_lst:
            view.addResLabels(reslabel[0], reslabel[1])


    if zoom_dict is None:
        view.zoomTo()
    else:
        view.zoomTo(zoom_dict)

    view.spin(spin_on)

    view.zoom(zoom)
    showmol(view, height=height, width=width)


def show_st_structure(df,
                    mut_resids=None, stick_resids=None,
                    label_muts=None, label_resids=False, label_ligs=False, label_prots=False, 
                    cartoon_style="oval",
                    cartoon_trans=1, surface_trans=0, mut_trans=0.5,
                    mut_color=None, sw1_color=None, sw2_color=None,
                    aa_scheme=False,
                    spin_on=False,
                    all_chains=False,
                    zoom_resids=None,
                    zoom=1.5,
                    width=400,
                    height=400,
                    show_legend=True,
                    legend_font_size="medium",
                    st_col=None):

    if sw1_color is None:
        sw1_color = loop_color_dict[sw1_name]
    if sw2_color is None:
        sw2_color = loop_color_dict[sw2_name]
    if mut_color is None:
        mut_color = red_hex

    stick_resid_lst = res_to_lst(stick_resids)
    mut_resid_lst = res_to_lst(mut_resids)

    style_lst = list()
    label_lst = list()
    reslabel_lst = list()
    surface_lst = list()

    opacity = 0
    if all_chains:
        opacity = 0.75
    
    pdb_code = df[pdb_code_col].iloc[0]
    chainid = df[chainid_col].iloc[0]

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
                    "opacity": opacity,
                }
            },
        ]
    )

    style_lst.append(
        [
            {
                "chain": chainid,
            },
            {
                "cartoon": {
                    "color": "white",
                    "style": cartoon_style,
                    "thickness": 0.2,
                    "opacity": cartoon_trans,
                }
            },
        ]
    )

    surface_lst = [
        [
            {"opacity": surface_trans, "color": "white"},
            {"chain": chainid, "hetflag": False},
        ]
    ]

    if stick_resids is not None:
        if aa_scheme:
            style_lst.append(
                [
                    {"chain": chainid,
                    "resi": stick_resid_lst,
                    "elem": "C"},
                    {"stick": {"colorscheme": "amino", "radius": 0.2}},
                ]
            ) 
            style_lst.append(
                    [
                        {"chain": chainid, "resi": stick_resid_lst},
                        {"stick": {"radius": 0.2}},
                    ]
                )
            surface_lst.append(
            [
                {"opacity": surface_trans,"colorscheme": "amino"},
                {"chain": chainid, "resi": stick_resid_lst, "hetflag": False},
            ]
            )
        else:
            for stick_resid in stick_resid_lst:

                resid_color = "white"
                if int(stick_resid) in sw1_resid_lst:
                    resid_color = sw1_color
                if int(stick_resid) in sw2_resid_lst:
                    resid_color = sw2_color

                style_lst.append(
                [
                    {
                        "chain": chainid,
                        "resi": stick_resid,
                        "elem": "C",
                    },
                    {"stick": {"color": resid_color, "radius": 0.2}},
                ]
            )

                style_lst.append(
                    [
                        {"chain": chainid, "resi": stick_resid},
                        {"stick": {"radius": 0.2}},
                    ]
                )  

        for stick_resid in stick_resid_lst:

            add_reslabel = label_resids
            if type(label_resids) == dict:
                if stick_resid in list(label_resids.keys()):
                    add_reslabel = label_resids[stick_resid]

            if add_reslabel:   
                reslabel_lst.append(
                        [
                            {"chain": chainid,
                            "resi": stick_resid},
                            {
                                "backgroundColor": "lightgray",
                                "fontColor": "black",
                                "backgroundOpacity": 0.5,
                            },
                        ]
                    )


    if mut_resids is not None:
        for mut_resid in mut_resid_lst:
            style_lst.append(
            [
                {
                    "chain": chainid,
                    "resi": mut_resid,
                    "elem": "C",
                },
                {"stick": {"color": mut_color, "radius": 0.2}},
            ]
        )

            style_lst.append(
                [
                    {"chain": chainid, "resi": mut_resid},
                    {"stick": {"radius": 0.2}},
                ]
            )     

            surface_lst.append([
                    
                        {"opacity": mut_trans, "color": mut_color},
                        {"chain": chainid, "resi": mut_resid},
                    
                ])

            add_reslabel = label_muts
            if type(label_muts) == dict:
                if mut_resid in list(label_muts.keys()):
                    add_reslabel = label_muts[mut_resid]

            if add_reslabel:
                reslabel_lst.append(
                        [
                            {"chain": chainid,
                            "resi": mut_resid},
                            {
                                "backgroundColor": "lightgray",
                                "fontColor": "black",
                                "backgroundOpacity": 0.5,
                            },
                        ]
                    )
            
    if not aa_scheme:
        for loop_name, loop_resids in loop_resid_dict.items():

            if loop_name == sw1_name:
                loop_color = sw1_color
            elif loop_name == sw2_name:
                loop_color = sw2_color

            surface_lst.append(
                [
                    {"opacity": surface_trans, "color": loop_color},
                    {"chain": chainid, "resi": loop_resids, "hetflag": False},
                ]
            )

            style_lst.append(
                [
                    {
                        "chain": chainid,
                        "resi": loop_resids,
                    },
                    {
                        "cartoon": {
                            "style": cartoon_style,
                            "color": loop_color,
                            "thickness": 0.2,
                            "opacity": cartoon_trans,
                        }
                    },
                ]
            )

    for lig_col in lig_col_lst:

        lig_lst = str_to_lst(df[lig_col].iloc[0])

        if "None" not in lig_lst:

            lig_style = "stick"
            lig_color = "whiteCarbon"
            lig_scheme = "colorscheme"
            lig_radius = 0.2

            if lig_col == ion_lig_col:

                lig_style = "sphere"
                lig_scheme = "color"
                lig_color = "chartreuse"
                lig_radius = 0.8

            for lig in lig_lst:

                if lig_col != pharm_lig_col:

                    lig_sele = {
                                "resn": lig,
                            }

                    if lig_col != mem_lig_col:
                         lig_sele["chain"] = chainid

                    style_lst.append(
                        [
                            lig_sele,
                            {
                                lig_style: {
                                    lig_scheme: lig_color,
                                    "radius": lig_radius,
                                }
                            },
                        ]
                    )

                add_reslabel = label_ligs
                if type(label_ligs) == dict:
                    if lig in list(label_ligs.keys()):
                        add_reslabel = label_ligs[lig]
                
                if add_reslabel:
                    reslabel_lst.append(
                        [
                            {
                                "chain": chainid,
                                "resn": lig,
                            },
                            {
                                "backgroundColor": "lightgray",
                                "fontColor": "black",
                                "backgroundOpacity": 0.5,
                            },
                        ]
                    )

    pharm_class_lst = str_to_lst(df[pharm_class_col].iloc[0])

    if none_pharm_name not in pharm_class_lst:
        if len(pharm_class_lst) > 1:
            pharm_class = other_pharm_name
        else:
            pharm_class = pharm_class_lst[0]

        pharm_color = pharm_color_dict[pharm_class]

        pharm_lig_lst = str_to_lst(df.at[0, pharm_lig_col])

        for pharm_lig in pharm_lig_lst:
            style_lst.append(
                [
                    {
                        "chain": chainid,
                        "resn": pharm_lig,
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
                    },
                    {"stick": {"radius": 0.2}},
                ]
            )

            if pharm_class == sp2_name:
                style_lst.append(
                    [
                        {
                            "chain": chainid,
                            "resi": 12,
                        },
                        {"stick": {"radius": 0.2}},
                    ]
                )

    prot_class_lst = str_to_lst(df[prot_class_col].iloc[0])

    if none_prot_name not in prot_class_lst:
        if len(prot_class_lst) > 1:
            prot_class = other_prot_name
        else:
            prot_class = prot_class_lst[0]
        prot_color = prot_color_dict[prot_class]

        bound_chainid_lst = str_to_lst(df[bound_prot_chainid_col].iloc[0])

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

            add_chainlabel = label_prots
            if type(label_prots) == dict:
                if bound_prot_chainid in list(label_prots.keys()):
                    add_chainlabel = label_prots[bound_prot_chainid]

            if add_chainlabel:

                prot_label = f"Chain {bound_prot_chainid}"

                label_lst.append(
                    [
                        prot_label,
                        {
                            "backgroundColor": "lightgray",
                            "fontColor": "black",
                            "backgroundOpacity": 0.8,
                        },
                        {"chain": bound_prot_chainid},
                    ]
                )

    zoom_dict = {"chain": chainid}

    if zoom_resids is not None:
        if type(zoom_resids) == dict:
            zoom_dict = merge_dicts([zoom_dict, zoom_resids])
        else:
            zoom_dict["resi"] = zoom_resids

    end_str = f"Powered by [Stmol](https://github.com/napoles-uach/streamlit_3dmol) (PDB: [{pdb_code.upper()}](https://www.rcsb.org/structure/{pdb_code}))."

    if st_col is None:
        show_st_3dmol(
            pdb_code,
            style_lst=style_lst,
            surface_lst=surface_lst,
            reslabel_lst=reslabel_lst,
            label_lst=label_lst,
            cartoon_style=cartoon_style,
            spin_on=spin_on,
            zoom_dict=zoom_dict,
            zoom=zoom,
            width=width,
            height=height,
        )
        st.markdown(end_str, unsafe_allow_html=True)

    else:
        with st_col:
            show_st_3dmol(
                pdb_code,
                style_lst=style_lst,
                surface_lst=surface_lst,
                reslabel_lst=reslabel_lst,
                label_lst=label_lst,
                cartoon_style=cartoon_style,
                spin_on=spin_on,
                zoom_dict=zoom_dict,
                zoom=zoom,
                width=width,
                height=height,
            )                 
        st_col.markdown(end_str, unsafe_allow_html=True)

    if show_legend:
        if aa_scheme:
            aa_type_dict = {
                    'Acidic':['ASP','GLU'],
                    'Basic':['LYS', 'ARG','HIS'],
                    'Polar':['ASN','GLN','SER','THR', 'CYS'],
                    'Nonpolar':['ILE', 'VAL', 'LEU', 'MET','PRO', 'GLY', 'ALA'],
                    'Aromatic':['PHE', 'TYR', 'TRP'],
                    }

            aa_color_dict = {'ASP':'#E60A0A','GLU':'#E60A0A','CYS':'#E6E600','MET':'#E6E600',
                            'LYS':'#145AFF','ARG':'#145AFF','SER':'#FA9600','THR':'#FA9600',
                            'PHE':'#3232AA','TYR':'#3232AA','ASN':'#00DCDC','GLN':'#00DCDC',
                            'GLY':'#C8C8C8','LEU':'#0F820F','VAL':'#0F820F','ILE':'#0F820F',
                            'ALA':'#C8C8C8','TRP':'#B45AB4','HIS':'#8282D2','PRO':'#DC9682'}

            for aa_type, aa_lst in aa_type_dict.items():

                legend_str = get_html_text({f"{aa_type}: ": "#31333F"}, font_weight='bold', font_size=legend_font_size)

                for i, aa_name in enumerate(aa_lst): 

                    legend_str += get_html_text({aa_name:aa_color_dict[aa_name]}, font_size=legend_font_size)

                    if i != len(aa_lst) - 1:
                        legend_str += get_html_text({", ":"#31333F"}, font_size=legend_font_size)

                if st_col is None:
                    st.markdown(legend_str, unsafe_allow_html=True)
                else:
                    st_col.markdown(legend_str, unsafe_allow_html=True)
       
            
            


     




    
