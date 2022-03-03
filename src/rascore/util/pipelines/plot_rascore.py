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

from ..scripts.build_dih_table import build_dih_table
from ..scripts.build_dist_table import build_dist_table
from ..scripts.make_facet_plot import make_facet_plot, grid_hex
from ..scripts.write_pymol_script import write_pymol_script
from ..scripts.prep_pocket import pocket_bound_name, pocket_unbound_name

from ..constants.dimer import dimer_name
from ..constants.prot import (
    prot_class_lst,
    prot_color_dict,
    none_prot_name,
    gap_name,
    gef_cdc_name,
    binder_name,
    nano_name,
    other_prot_name,
    mult_prot_name
)
from ..constants.pharm import (
    sp12_name,
    sp2_name,
    other_pharm_name,
    pharm_color_dict,
    pharm_class_lst,
    none_pharm_name,
    other_pharm_name,
    mult_pharm_name
)
from ..constants.gene import hras_name
from ..constants.nuc import nuc_class_lst, nuc_color_dict, gtp_atomids
from ..constants.conf import (
    loop_resid_dict,
    conf_color_dict,
    conf_nuc_color_dict,
    outlier_name,
    disorder_name,
    sw1_name,
    sw2_name,
    sw1_resids,
    sw2_resids,
    sw1_color,
    sw2_color,
    y32_name,
    q61_name,
    y71_name,
    resid_color_dict,
    sw1_gdp_out_off_name,
    sw1_gtp_in_on_name,
    sw1_nf_out_gef_name,
    sw2_gtp_in_r_name,
    sw2_gtp_in_sp12a_name,
    sw2_gdp_in_sp12_name,
    sw2_gdp_out_binder_name,
    sw2_gdp_out_sp2a_name,
    sw2_gdp_out_sp2b_name,
    sw2_gtp_in_sp12b_name,
    sw2_gtp_out_t_name,
    sw2_nf_out_gef_name
)
from ..constants.pocket import pocket_site_lst, pocket_color_dict
from ..constants.pml import (
    sup_resids,
    show_resids,
    sup_pdb_code,
    sup_chainid,
    mono_view,
    y32_view,
    gap_view,
    gef_view,
    mut_view,
    sw1_mono_view,
    sw2_mono_view,
    sw1_prot_view,
    sw2_prot_view,
)
from ..functions.stat import calc_rr
from ..functions.plot import make_venn_plot, make_stacked_barplot
from ..functions.color import change_hex_alpha, get_palette_hex_lst, gray_hex, blue_hex, orange_hex, green_hex, purple_hex, pink_hex, cyan_hex
from ..functions.file import (
    entry_table_file,
    dist_table_file,
    dih_table_file,
    pocket_table_file,
    interf_table_file,
    dih_json_file,
    pymol_pml_file,
    plot_img_file,
    venn_img_file,
    stat_table_file
)
from ..functions.path import (
    get_file_path,
    load_table,
    load_json,
    save_table,
    get_core_path,
    rascore_str,
    plot_str,
)
from ..functions.col import (
    rama_col,
    nuc_class_col,
    phi_col,
    psi_col,
    cluster_col,
    pdb_id_col,
    core_path_col,
    bio_lig_col,
    atom_dist_col,
    hb_status_col,
    outlier_col,
    gene_class_col,
    mut_class_col,
    chi1_col,
    chi2_col,
    pocket_class_col,
    pocket_type_col,
    pocket_site_col,
    pocket_status_col,
    pocket_volume_col,
    pocket_score_col,
    pocket_path_col,
    pharm_class_col,
    rotamer_col,
    match_class_col,
    prot_class_col,
    interf_path_col,
    bound_interf_chainid_col,
    bound_prot_cont_col,
    bound_lig_cont_col,
    sig_col,
    risk_ratio_col
)

from ..functions.table import (
    mask_equal,
    mask_greater,
    mask_unequal,
    get_col_most_common,
    lst_col,
    build_col_count_dict,
    make_dict,
)
from ..functions.lst import res_to_lst, str_to_lst


def plot_rama(df, dih_dict, out_path, num_cpu=1):

    dih_table_path = get_file_path(
        dih_table_file, dir_str=rama_col, dir_path=out_path, pre_str=False
    )

    g_resids = "2-166"

    # build_dih_table(
    #     df, dih_dict, dih_table_path=dih_table_path, bb_resids=g_resids, num_cpu=num_cpu
    # )

    # dih_df = load_table(dih_table_path)

    # make_facet_plot(
    #     dih_df,
    #     get_file_path(f"{rama_col}_{plot_img_file}", dir_str=rama_col, dir_path=out_path, pre_str=False),
    #     x_col=phi_col,
    #     y_col=psi_col,
    #     x_str="φ",
    #     y_str="ψ",
    #     hue_col=nuc_class_col,
    #     hue_palette=nuc_color_dict,
    #     col_wrap=11,
    #     plot_width=5,
    #     plot_height=7,
    #     darken_lst=[
    #         res_to_lst(sup_resids),
    #     ],
    #     darken_palette=[grid_hex],
    #     rename_col=res_to_lst(g_resids),
    #     show_legend=True,
    #     hue_count=True,
    #     legend_pad=5,
    # )

    major_name = "Major"
    minor_name = "Minor"

    color_dict = dict()

    for nuc_class in nuc_class_lst:
        nuc_color = nuc_color_dict[nuc_class]
        color_dict[f"{nuc_class}-{major_name}"] = nuc_color
        color_dict[f"{nuc_class}-{minor_name}"] = change_hex_alpha(nuc_color, 0.5)

    for loop_name, loop_resids in loop_resid_dict.items():

        dih_table_path = get_file_path(
            dih_table_file, dir_str=loop_name, dir_path=out_path
        )

        # build_dih_table(
        #     df,
        #     dih_dict,
        #     dih_table_path=dih_table_path,
        #     bb_resids=loop_resids,
        #     num_cpu=num_cpu,
        # )

        dih_df = load_table(dih_table_path)

        rama_dict = dict()

        for cluster in lst_col(dih_df, loop_name, unique=True):
            if cluster not in [outlier_name, disorder_name]:
                rama_dict[cluster] = [
                    x
                    for x in get_col_most_common(
                        mask_equal(dih_df, loop_name, cluster), rama_col
                    )
                    if "-" not in x
                ][0]

        for index in list(dih_df.index.values):
            cluster = dih_df.at[index, loop_name]
            rama_status = major_name
            if cluster not in [outlier_name, disorder_name]:
                if dih_df.at[index, rama_col] != rama_dict[cluster]:
                    rama_status = minor_name
            dih_df.at[index, nuc_class_col] += "-"
            dih_df.at[index, nuc_class_col] += rama_status

        make_facet_plot(
            dih_df,
            get_file_path(
                f"{rama_col}_{plot_img_file}", dir_str=loop_name, dir_path=out_path
            ),
            x_col=phi_col,
            y_col=psi_col,
            x_str="φ",
            y_str="ψ",
            hue_col=nuc_class_col,
            row_col=loop_name,
            row_order=list(conf_nuc_color_dict[loop_name].keys()),
            rename_col=res_to_lst(loop_resids),
            hue_palette=color_dict,
            row_palette=conf_nuc_color_dict[loop_name],
            plot_width=5,
            show_legend=False,
            marker_size=1,
        )


def plot_pymol(df, interf_df, sup_core_path, out_path):

    dark_hex_lst = get_palette_hex_lst('Dark2')

    sw1_pymol_name_dict = {'1q21A':"GDP-bound","1bkdR":"Nucleotide-Free","4eflA":"State 1","5p21A":"State 2"}
    sw2_pymol_name_dict = {"5b2zA":"T state","3k8yA":"R State"}

    sw1_pymol_color_dict = {"GDP-bound":dark_hex_lst[0],"Nucleotide-Free":dark_hex_lst[1],"State 1":dark_hex_lst[2],"State 2":dark_hex_lst[3]}
    sw2_pymol_color_dict = {"T state":dark_hex_lst[4],"R State":dark_hex_lst[5]}

    sw1_df = mask_equal(df, pdb_id_col, list(sw1_pymol_name_dict.keys()))
    sw2_df = mask_equal(df, pdb_id_col, list(sw2_pymol_name_dict.keys()))

    sw1_df[sw1_name] = sw1_df[pdb_id_col].map(sw1_pymol_name_dict)
    sw2_df[sw2_name] = sw2_df[pdb_id_col].map(sw2_pymol_name_dict)

    write_pymol_script(
        sw1_df,
        get_file_path(f"{sw1_name}_{pymol_pml_file}", dir_path=out_path),
        stick_resids=[32],
        loop_resids=[sw1_resids],
        group_col=sw1_name,
        color_palette=sw1_pymol_color_dict,
        color_group=True,
        show_bio=True,
        show_ion=True,
        style_ribbon=False,
        sup_resids=sup_resids,
        sup_coord_path=sup_core_path,
        sup_chainid=sup_chainid,
        set_view=mono_view,
        show_resids=show_resids,
    )

    write_pymol_script(
        sw2_df,
        get_file_path(f"{sw2_name}_{pymol_pml_file}", dir_path=out_path),
        stick_resids=[71],
        loop_resids=[sw2_resids],
        group_col=sw2_name,
        color_palette=sw2_pymol_color_dict,
        color_group=True,
        show_bio=True,
        show_ion=True,
        style_ribbon=False,
        sup_resids=sup_resids,
        sup_coord_path=sup_core_path,
        sup_chainid=sup_chainid,
        set_view=mono_view,
        show_resids=show_resids,
    )

    for resid_name in [y32_name, y71_name]:

        if resid_name == y32_name:
            loop_resids = sw1_resids
            stick_resids = [32]
        elif resid_name == y71_name:
            loop_resids = sw2_resids
            stick_resids = [71]

        for nuc_class in nuc_class_lst:

            write_pymol_script(
                mask_equal(mask_unequal(df,resid_name,disorder_name),nuc_class_col,nuc_class),
                get_file_path(f"{nuc_class}_{pymol_pml_file}", dir_str=resid_name, dir_path=out_path),
                group_col=resid_name,
                stick_resids=stick_resids,
                style_ribbon=True,
                thick_bb=False,
                show_bio=True,
                color_palette=resid_color_dict[resid_name],
                sup_resids=sup_resids,
                sup_coord_path=sup_core_path,
                sup_chainid=sup_chainid,
                set_view=mono_view,
                show_resids=show_resids,
            )


    for loop_name, loop_resids in loop_resid_dict.items():

        if loop_name == sw1_name:
            stick_resids = [32]
        elif loop_name == sw2_name:
            stick_resids = [71]

        write_pymol_script(
            mask_unequal(df,loop_name,[outlier_name,disorder_name]),
            get_file_path(pymol_pml_file, dir_str=loop_name, dir_path=out_path),
            group_col=loop_name,
            stick_resids=stick_resids,
            loop_resids=loop_resids,
            style_ribbon=True,
            thick_bb=False,
            show_bio=True,
            color_palette=conf_color_dict[loop_name],
            sup_resids=sup_resids,
            sup_coord_path=sup_core_path,
            sup_chainid=sup_chainid,
            set_view=mono_view,
            show_resids=show_resids,
        )

    for loop_name, loop_resids in loop_resid_dict.items():
        for prot_class in [x for x in prot_class_lst if x != none_prot_name]:

            prot_df = mask_equal(df, prot_class_col, prot_class)

            if loop_name == sw1_name:
                stick_resids = [32]
                if prot_class == gap_name:
                    set_view = gap_view
                elif gef_cdc_name in prot_class:
                    set_view = gef_view
                elif binder_name in prot_class:
                    set_view = sw1_mono_view
                else:
                    set_view = sw1_prot_view
            elif loop_name == sw2_name:
                stick_resids = [71]
                set_view = sw2_prot_view
                if binder_name in prot_class:
                    set_view = sw2_mono_view

            write_pymol_script(
                prot_df,
                get_file_path(
                    f"{loop_name}_{prot_class}_{pymol_pml_file}",
                    dir_str=prot_class_col,
                    dir_path=out_path,
                    pre_str=False
                ),
                group_col=loop_name,
                stick_resids=stick_resids,
                loop_resids=loop_resids,
                style_ribbon=True,
                thick_bb=False,
                show_bio=True,
                show_prot=prot_color_dict[prot_class.split(".")[0]],
                color_palette=conf_color_dict[loop_name],
                sup_resids=sup_resids,
                sup_coord_path=sup_core_path,
                sup_chainid=sup_chainid,
                set_view=set_view,
                show_resids=show_resids,
            )

    for loop_name, loop_resids in loop_resid_dict.items():
        for pharm_class in [x for x in pharm_class_lst if x != none_pharm_name]:

            pharm_df = mask_equal(df, pharm_class_col, pharm_class)

            if loop_name == sw1_name:
                stick_resids = [32]
                set_view = sw1_mono_view
            elif loop_name == sw2_name:
                stick_resids = [71]
                set_view = sw2_mono_view

            if sp2_name in pharm_class:
                stick_resids += [12]

            write_pymol_script(
                pharm_df,
                get_file_path(
                    f"{loop_name}_{pharm_class}_{pymol_pml_file}",
                    dir_str=pharm_class_col,
                    dir_path=out_path,
                ),
                group_col=loop_name,
                stick_resids=stick_resids,
                loop_resids=loop_resids,
                thick_bb=False,
                show_bio=True,
                show_pharm=pharm_color_dict[pharm_class],
                color_palette=conf_color_dict[loop_name],
                sup_resids=sup_resids,
                sup_coord_path=sup_core_path,
                sup_chainid=sup_chainid,
                set_view=set_view,
                show_resids=show_resids,
            )

    for loop_name, loop_resids in loop_resid_dict.items():

        write_pymol_script(
            interf_df,
            get_file_path(
                f"{loop_name}_{pymol_pml_file}", dir_str=dimer_name, dir_path=out_path
            ),
            loop_resids=loop_resids,
            group_col=loop_name,
            color_palette=conf_color_dict[loop_name],
            style_ribbon=True,
            thick_bb=False,
            show_bio=True,
            show_prot=True,
            coord_path_col=interf_path_col,
            prot_chainid_col=bound_interf_chainid_col,
            sup_resids=sup_resids,
            sup_coord_path=sup_core_path,
            sup_chainid=sup_chainid,
            set_view=mono_view,
            show_resids=show_resids,
        )


def plot_dist(df, dih_dict, sup_core_path, out_path):

    dist_table_path = get_file_path(
        f"{y32_name}_{y71_name}_{dist_table_file}", dir_path=out_path
    )

    # build_dist_table(
    #         df,
    #         x_resids=[32, 71],
    #         y_resids=[12, 9],
    #         x_atomids=['OH','OH'],
    #         y_atomids=['CA','CA'],
    #         atom_dist_col_lst=[y32_name, y71_name],
    #         dist_table_path=dist_table_path
    #     )

    dist_df = load_table(dist_table_path)

    make_facet_plot(
        dist_df,
        get_file_path(f"{y32_name}_{plot_img_file}", dir_str=y32_name, dir_path=out_path,pre_str=False),
        x_col=y32_name,
        x_str="Y32(OH):G12(CA)\nDistance (Å)",
        plot_width=1.5,
        plot_height=1.5,
        plot_kde=True,
        kde_bw=0.3,
        show_legend=False,
        v_lines=[10.5],
        x_round=0,
        y_round=1,
        x_lim=[0,20],
        x_ticks=[0,5,10,15,20],
        hue_col=nuc_class_col,
        hue_palette=nuc_color_dict
    )

    make_facet_plot(
        dist_df,
        get_file_path(f"{y71_name}_{plot_img_file}", dir_str=y71_name, dir_path=out_path,pre_str=False),
        x_col=y71_name,
        x_str="Y71(OH):V9(CA)\nDistance (Å)",
        plot_width=1.5,
        plot_height=1.5,
        plot_kde=True,
        kde_bw=0.3,
        show_legend=False,
        v_lines=[8.75],
        x_round=0,
        y_round=1,
        x_lim=[0,20],
        x_ticks=[0,5,10,15,20],
        hue_col=nuc_class_col,
        hue_palette=nuc_color_dict
    )

    # sw1_gtp_df = mask_unequal(df, sw1_name, [sw1_nf_name, sw1_gdp_name, outlier_name, disorder_name])

    # dist_table_path = get_file_path(
    #     dist_table_file, dir_str=y32_name, dir_path=out_path
    # )

    # build_dist_table(
    #     sw1_gtp_df,
    #     x_resids=[32],
    #     y_resids=[bio_lig_col],
    #     x_atomids=["OH"],
    #     y_atomids=[gtp_atomids],
    #     atom_dist_col_lst=[atom_dist_col],
    #     hb_status_col_lst=[hb_status_col],
    #     outlier_col_lst=[outlier_col],
    #     dist_table_path=dist_table_path,
    #     check_hb=True,
    # )

    # dist_df = load_table(dist_table_path)

    # dist_df[hb_status_col] = dist_df[hb_status_col].map(sw1_gtp_dict)

    # make_facet_plot(
    #     dist_df,
    #     get_file_path(plot_img_file, dir_str=y32_name, dir_path=out_path),
    #     x_col=atom_dist_col,
    #     x_str="Y32(OH):3P(O1G)\nDistance (Å)",
    #     plot_width=1.5,
    #     plot_height=1.5,
    #     plot_kde=True,
    #     show_legend=False,
    #     x_round=1,
    #     y_round=1,
    #     x_ticks=[2.5, 5.0, 7.5, 10.0],
    #     hue_col=hb_status_col,
    #     hue_palette=sw1_gtp_color_dict,
    # )

    # hras_df = mask_equal(
    #     mask_equal(mask_equal(dist_df, bio_lig_col, "GNP"), outlier_col, "False"),
    #     gene_class_col,
    #     hras_name,
    # )
    # write_pymol_script(
    #     hras_df,
    #     get_file_path(
    #         f"{hras_name}_{pymol_pml_file}", dir_str=y32_name, dir_path=out_path
    #     ),
    #     stick_resids=[32],
    #     group_col=hb_status_col,
    #     color_palette=sw1_gtp_color_dict,
    #     thick_bb=False,
    #     show_bio=True,
    #     sup_coord_path=sup_core_path,
    #     sup_chainid=sup_chainid,
    #     set_view=y32_view,
    #     x_hb_resids=[32],
    #     x_hb_atomids=["OH"],
    #     y_hb_resids=[bio_lig_col],
    #     y_hb_atomids=[gtp_atomids],
    #     show_hb=[sw1_gtp_dir_color],
    #     show_wmhb=[sw1_gtp_wat_color],
    #     show_resids=show_resids,
    # )

    # r_df = mask_equal(dist_df, sw2_name, sw2_gtp_r_name)

    # dih_table_path = get_file_path(
    #     f"{q61_name}_{dih_table_file}",
    #     dir_str=y32_name,
    #     dir_path=out_path,
    # )

    # build_dih_table(
    #     r_df, dih_dict, chi1_resids=61, chi2_resids=61, dih_table_path=dih_table_path
    # )

    # dih_df = load_table(dih_table_path)

    # for mut_class in ["WT", "G12D", "G12V"]:

    #     mut_df = mask_equal(
    #         dih_df,
    #         mut_class_col,
    #         mut_class,
    #     )

    #     make_facet_plot(
    #         mut_df,
    #         get_file_path(
    #             f"{mut_class}_{q61_name}_{chi1_col}_{chi2_col}_{plot_img_file}",
    #             dir_str=y32_name,
    #             dir_path=out_path,
    #         ),
    #         x_col=chi1_col,
    #         y_col=chi2_col,
    #         hue_col=hb_status_col,
    #         hue_palette=sw1_gtp_color_dict,
    #         x_str=f"χ$^1$",
    #         y_str=f"χ$^2$",
    #         plot_width=1.5,
    #         plot_height=1.5,
    #         kde_bw=1,
    #         show_legend=False,
    #     )

    #     write_pymol_script(
    #         mut_df,
    #         get_file_path(
    #             f"{mut_class}_{pymol_pml_file}", dir_str=y32_name, dir_path=out_path
    #         ),
    #         stick_resids=[12, 32, 61],
    #         group_col=hb_status_col,
    #         color_palette=sw1_gtp_color_dict,
    #         thick_bb=False,
    #         show_bio=True,
    #         sup_coord_path=sup_core_path,
    #         sup_chainid=sup_chainid,
    #         set_view=mut_view,
    #         x_hb_resids=[12, 12, 32, 61],
    #         x_hb_atomids=["OD1", "OD1", "OH", "OE1"],
    #         y_hb_resids=[32, bio_lig_col, bio_lig_col, 32],
    #         y_hb_atomids=["OH", gtp_atomids, gtp_atomids, "OH"],
    #         show_hb=[
    #             sw1_gtp_wat_color,
    #             sw1_gtp_wat_color,
    #             sw1_gtp_dir_color,
    #             sw1_gtp_dir_color,
    #         ],
    #         show_wmhb=[sw1_gtp_wat_color] * 4,
    #         show_resids=show_resids,
    #     )


def plot_pockets(df, pocket_df, sup_core_path, out_path):

    make_facet_plot(
        pocket_df,
        get_file_path(
            f"{pocket_volume_col}_{pocket_score_col}_{plot_img_file}",
            dir_str=pocket_class_col,
            dir_path=out_path,
        ),
        x_col=pocket_volume_col,
        y_col=pocket_score_col,
        x_str="Pocket Volume (Å$^3$)",
        y_str="Druggability Score",
        y_round=1,
        col_col=pocket_site_col,
        col_order=pocket_site_lst,
        col_wrap=3,
        hue_col=pocket_type_col,
        hue_palette=pocket_color_dict,
        hue_order=list(pocket_color_dict.keys()),
        show_legend=False,
        plot_width=2.6,
        plot_height=1.2,
        x_rotation=45,
        x_ha="right",
        plot_reg=True,
        log_reg=True,
        trun_reg=True,
        x_lim=[0, 1700],
        x_ticks=[0, 500, 1000, 1500],
        y_lim=[0, 1.2],
        y_ticks=[0.0, 0.5, 1.0],
    )

    for pocket_site in pocket_site_lst:

        if pocket_site == other_pharm_name:
            show_pocket = True
        else:
            show_pocket = pocket_color_dict[pocket_site]

        stick_resids = None
        if pocket_site == sp2_name:
            stick_resids = [12]

        fp_df = mask_equal(pocket_df, pocket_site_col, pocket_site)

        write_pymol_script(
            fp_df,
            get_file_path(
                f"{pocket_site}_{pymol_pml_file}",
                dir_str=pocket_class_col,
                dir_path=out_path,
            ),
            stick_resids=stick_resids,
            loop_resids=[sw1_resids, sw2_resids],
            color_palette=[sw1_color, sw2_color],
            group_col=pocket_type_col,
            color_group=False,
            thick_bb=False,
            sup_resids=sup_resids,
            sup_coord_path=sup_core_path,
            sup_chainid=sup_chainid,
            show_bio=True,
            show_pocket=show_pocket,
            set_view=mono_view,
            coord_path_col=pocket_path_col,
            show_resids=show_resids,
        )

        if pocket_site != other_pharm_name:
            all_df = mask_equal(df, pharm_class_col, pocket_site)

            del fp_df[pocket_path_col]

            fpocket_lst = lst_col(fp_df, pdb_id_col, unique=True)
            all_lst = lst_col(all_df, pdb_id_col, unique=True)

            make_venn_plot(
                fpocket_lst,
                all_lst,
                get_file_path(
                    f"{pocket_site}_{venn_img_file}",
                    dir_str=pocket_class_col,
                    dir_path=out_path,
                ),
                color_1=gray_hex,
                color_2=pocket_color_dict[pocket_site],
                label_1="Predicted",
                label_2="Observed",
                plot_height=1.5,
                plot_width=1.5,
                alpha=0.75,
            )

            write_pymol_script(
                all_df,
                get_file_path(
                    f"{pocket_site}_{y71_name}_{pymol_pml_file}",
                    dir_str=pocket_class_col,
                    dir_path=out_path,
                ),
                stick_resids=[71],
                loop_resids=[71],
                color_palette=resid_color_dict[y71_name],
                group_col=y71_name,
                thick_bb=False,
                sup_resids=sup_resids,
                sup_coord_path=sup_core_path,
                sup_chainid=sup_chainid,
                show_bio=True,
                show_pharm=show_pocket,
                set_view=mono_view,
                show_resids=show_resids,
            )

    pocket_loop_cluster_dict = {
        sp12_name: {
            sw2_name: [
                sw2_gtp_in_r_name,
                sw2_gtp_in_sp12a_name,
                sw2_gtp_in_sp12b_name,
                sw2_gdp_in_sp12_name,
                outlier_name,
                disorder_name
            ]
        },
        sp2_name: {
            sw2_name: [
                sw2_gdp_out_sp2a_name,
                sw2_gdp_out_sp2b_name,
                outlier_name,
                disorder_name
            ]
        },
    }

    pocket_loop_cluster_dict[sp12_name][sw1_name] = [
        sw1_gtp_in_on_name,
        sw1_gdp_out_off_name,
        outlier_name,
        disorder_name
    ]
    pocket_loop_cluster_dict[sp2_name][sw1_name] = [
        sw1_gdp_out_off_name,
        outlier_name,
        disorder_name
    ]

    for pocket_site in [sp2_name, sp12_name]:
        all_df = mask_equal(df, pharm_class_col, pocket_site)
        for loop_name in list(pocket_loop_cluster_dict[pocket_site].keys()):
            color_dict = conf_color_dict[loop_name].copy()

            if loop_name == sw1_name:
                stick_resids = [32]
            elif loop_name == sw2_name:
                stick_resids = [71]

            if pocket_site == sp2_name:
                stick_resids.append(12)

            write_pymol_script(
                all_df,
                get_file_path(
                    f"{loop_name}_{pocket_site}_{pymol_pml_file}",
                    dir_str=pocket_class_col,
                    dir_path=out_path,
                ),
                stick_resids=stick_resids,
                loop_resids=[loop_resid_dict[loop_name]],
                color_palette=color_dict,
                group_col=loop_name,
                color_group=True,
                thick_bb=False,
                sup_resids=sup_resids,
                sup_coord_path=sup_core_path,
                sup_chainid=sup_chainid,
                show_bio=True,
                show_pharm=pocket_color_dict[pocket_site],
                set_view=mono_view,
                coord_path_col=core_path_col,
                show_resids=show_resids,
            )

            make_stacked_barplot(
                plot_df=all_df,
                col_col=match_class_col,
                hue_col=loop_name,
                plot_path=get_file_path(
                    f"{loop_name}_{pocket_site}_{match_class_col}_{plot_img_file}",
                    dir_str=pocket_class_col,
                    dir_path=out_path,
                ),
                hue_palette=color_dict,
                x_str=f"% Structures ({pocket_site} Bound)",
                y_str=f"Inhibitor Chemistry",
                plot_height=1.0,
                plot_width=1.5,
                show_legend=False,
                col_count=True,
                show_barh=True,
            )

            for pocket_status in [pocket_bound_name, pocket_unbound_name]:

                fp_df = mask_equal(
                    mask_equal(
                        pocket_df,
                        pocket_site_col,
                        pocket_site,
                    ),
                    pocket_status_col,
                    pocket_status,
                )

                count_dict = build_col_count_dict(fp_df, loop_name)

                cluster_order = [
                    x
                    for x in list(color_dict.keys())
                    if x in list(count_dict.keys()) and count_dict[x] > 3
                ]

                make_facet_plot(
                    fp_df,
                    get_file_path(
                        f"{loop_name}_{pocket_site}_{pocket_status}_{pocket_volume_col}_{plot_img_file}",
                        dir_str=pocket_class_col,
                        dir_path=out_path,
                        pre_str=False,
                    ),
                    x_col=loop_name,
                    x_order=cluster_order,
                    x_palette=color_dict,
                    x_count=False,
                    y_col=pocket_volume_col,
                    x_str=f"{loop_name} Conformation ({pocket_site}-{pocket_status})",
                    y_str="Pocket Volume (Å$^3$)",
                    show_legend=False,
                    plot_width=1.0,
                    plot_height=1.0,
                    plot_kind="box",
                    marker_size=2,
                    x_pad=35,
                    y_lim=[0, 1700],
                    y_ticks=[0, 500, 1000, 1500],
                )

                make_facet_plot(
                    fp_df,
                    get_file_path(
                        f"{loop_name}_{pocket_site}_{pocket_status}_{pocket_score_col}_{plot_img_file}",
                        dir_str=pocket_class_col,
                        dir_path=out_path,
                        pre_str=False,
                    ),
                    x_col=loop_name,
                    x_order=cluster_order,
                    x_palette=color_dict,
                    x_count=False,
                    y_col=pocket_score_col,
                    x_str=f"{loop_name} Conformation ({pocket_site}-{pocket_status})",
                    y_str="Druggability Score",
                    show_legend=False,
                    plot_width=1.0,
                    plot_height=1.0,
                    plot_kind="box",
                    marker_size=2,
                    y_round=1,
                    x_pad=35,
                    y_lim=[0, 1.2],
                    y_ticks=[0.0, 0.5, 1.0],
                )

# def plot_table(df, out_path):

#     for group_col in [prot_class_col, match_class_col]:

#         if group_col == prot_class_col:
#             cont_col = bound_prot_cont_col
#             mask_lst = [mult_prot_name, other_prot_name, none_prot_name, binder_name, nano_name]
#         elif group_col == match_class_col:
#             cont_col = bound_lig_cont_col
#             mask_lst = [mult_pharm_name, other_pharm_name, none_prot_name]

#         cont_df = mask_unequal(df,group_col,mask_lst)

#         stat_df = pd.DataFrame()

#         i = 0
        
#         for index in list(cont_df.index.values):

#             group = cont_df.at[index, group_col]

#             obj_lst = str_to_lst(cont_df.at[index,cont_col],sep_txt=";")

#             for obj in obj_lst:
#                 cont_lst = str_to_lst(obj.split(':')[1])
#                 for cont in cont_lst:
#                     stat_df.at[i , group_col] = group
#                     stat_df.at[i , cont_col] = cont
#                     i += 1

#         stat_df = calc_rr(stat_df, group_col, cont_col)

#         stat_df = mask_unequal(stat_df,sig_col,'ns')
#         stat_df = mask_greater(stat_df,risk_ratio_col,2)


#         stat_table_path = get_file_path(
#                         stat_table_file,
#                         dir_str=group_col,
#                         dir_path=out_path,
#                     )

#         save_table(stat_table_path, stat_df)

def plot_rascore(build_path, out_path=None, num_cpu=1):

    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{plot_str}"

    entry_table_path = get_file_path(entry_table_file, dir_path=build_path)
    pocket_table_path = get_file_path(pocket_table_file, dir_path=build_path)
    interf_table_path = get_file_path(interf_table_file, dir_path=build_path)
    dih_json_path = get_file_path(dih_json_file, dir_path=build_path)

    sup_core_path = get_core_path(sup_pdb_code, sup_chainid, dir_path=build_path)

    df = load_table(entry_table_path)
    pocket_df = load_table(pocket_table_path)
    interf_df = load_table(interf_table_path)
    dih_dict = load_json(dih_json_path)

    sw1_dict = make_dict(lst_col(df, pdb_id_col), lst_col(df, sw1_name))
    sw2_dict = make_dict(lst_col(df, pdb_id_col), lst_col(df, sw2_name))

    pocket_df[sw1_name] = pocket_df[pdb_id_col].map(sw1_dict)
    pocket_df[sw2_name] = pocket_df[pdb_id_col].map(sw2_dict)

    interf_df[sw1_name] = interf_df[pdb_id_col].map(sw1_dict)
    interf_df[sw2_name] = interf_df[pdb_id_col].map(sw2_dict)

    #plot_pymol(df, interf_df, sup_core_path, out_path)
    # plot_dist(df, dih_dict, sup_core_path, out_path)
    plot_pockets(df, pocket_df, sup_core_path, out_path)
    #plot_rama(df, dih_dict, out_path,num_cpu=num_cpu)
