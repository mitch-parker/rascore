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
)
from ..constants.pharm import (
    sp12_name,
    sp2_name,
    other_pharm_name,
    pharm_color_dict,
    pharm_class_lst,
    none_pharm_name,
)
from ..constants.gene import hras_name
from ..constants.nuc import nuc_class_lst, nuc_color_dict, gtp_atomids
from ..constants.conf import (
    loop_resid_dict,
    conf_color_dict,
    conf_nuc_color_dict,
    sw1_gtp_name,
    sw1_nf_name,
    sw1_gdp_name,
    noise_name,
    sw1_gtp_dict,
    sw1_gtp_color_dict,
    sw1_name,
    sw2_name,
    sw1_resids,
    sw2_resids,
    sw1_color,
    sw2_color,
    sw1_gtp_dir_color,
    sw1_gtp_wat_color,
    sw1_gtp_wat_name,
    sw2_gtp_r_name,
    sw2_gtp_sp12_a_name,
    sw2_gtp_sp12_b_name,
    sw2_gdp_sp12_name,
    sw2_gdp_sp2_a_name,
    sw2_gdp_sp2_b_name,
    y32_name,
    q61_name,
    y71_name,
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
from ..functions.plot import make_venn_plot, make_stacked_barplot
from ..functions.color import change_hex_alpha, gray_hex
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
)
from ..functions.path import (
    get_file_path,
    load_table,
    load_json,
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
)

from ..functions.table import (
    mask_equal,
    mask_unequal,
    get_col_most_common,
    lst_col,
    build_col_count_dict,
    make_dict,
)
from ..functions.lst import res_to_lst


def plot_rama(df, dih_dict, out_path, num_cpu=1):

    dih_table_path = get_file_path(
        dih_table_file, dir_str=rama_col, dir_path=out_path, pre_str=False
    )

    g_resids = "2-166"

    build_dih_table(
        df, dih_dict, dih_table_path=dih_table_path, bb_resids=g_resids, num_cpu=num_cpu
    )

    dih_df = load_table(dih_table_path)

    make_facet_plot(
        dih_df,
        get_file_path(f"{rama_col}_{plot_img_file}", dir_path=out_path),
        x_col=phi_col,
        y_col=psi_col,
        x_str="φ",
        y_str="ψ",
        hue_col=nuc_class_col,
        hue_palette=nuc_color_dict,
        col_wrap=11,
        plot_width=5,
        plot_height=7,
        darken_lst=[
            res_to_lst(sup_resids),
        ],
        darken_palette=[grid_hex],
        rename_col=res_to_lst(g_resids),
        show_legend=True,
        hue_count=True,
        legend_pad=5,
    )

    major_name = "Major"
    minor_name = "Minor"

    color_dict = dict()

    for nuc_class in nuc_class_lst:
        nuc_color = nuc_color_dict[nuc_class]
        color_dict[f"{nuc_class}_{major_name}"] = nuc_color
        color_dict[f"{nuc_class}_{minor_name}"] = change_hex_alpha(nuc_color, 0.5)

    for loop_name, loop_resids in loop_resid_dict.items():

        dih_table_path = get_file_path(
            dih_table_file, dir_str=loop_name, dir_path=out_path
        )

        build_dih_table(
            df,
            dih_dict,
            dih_table_path=dih_table_path,
            bb_resids=loop_resids,
            num_cpu=num_cpu,
        )

        dih_df = load_table(dih_table_path)

        for index in list(dih_df.index.values):
            cluster = dih_df.at[index, loop_name]
            if sw1_gtp_name in cluster:
                cluster = sw1_gtp_name
            dih_df.at[index, loop_name] = cluster

        rama_dict = dict()

        for cluster in lst_col(dih_df, loop_name, unique=True):
            common_rama = [
                x
                for x in get_col_most_common(
                    mask_equal(dih_df, loop_name, cluster), rama_col
                )
                if "-" not in x
            ]
            rama_dict[common_rama]

        for index in list(dih_df.index.values):
            cluster = dih_df.at[index, loop_name]
            rama_status = major_name
            if cluster != "Noise":
                if dih_df.at[index, rama_col] != rama_dict[cluster]:
                    rama_status = minor_name
            dih_df.at[index, nuc_class_col] += "_"
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
            plot_width=5.5,
            show_legend=False,
            marker_size=1,
        )


def plot_pymol(df, interf_df, sup_core_path, out_path):

    write_pymol_script(
        mask_equal(df, core_path_col, sup_core_path),
        get_file_path(f"{sw1_name}_{sw2_name}_{pymol_pml_file}", dir_path=out_path),
        stick_resids=[32, 71],
        loop_resids=[sw1_resids, sw2_resids],
        group_col=nuc_class_col,
        color_palette=[sw1_color, sw2_color],
        color_group=False,
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
        df,
        get_file_path(f"{nuc_class_col}_{pymol_pml_file}", dir_path=out_path),
        loop_resids=[sw1_resids, sw2_resids],
        show_bio=True,
        style_ribbon=True,
        thick_bb=False,
        group_col=nuc_class_col,
        color_palette=[sw1_color, sw2_color],
        color_group=False,
        sup_resids=sup_resids,
        sup_coord_path=sup_core_path,
        sup_chainid=sup_chainid,
        set_view=mono_view,
        show_resids=show_resids,
    )

    gdp_b_lst = ["6bofA", "6bofB", "6m9wA", "6mqgA"]
    gtp_1_lst = ["1xcmA", "3kknA", "4eflA", "4efmA", ",4efnA", "6bp1A"]

    write_pymol_script(
        mask_equal(df, pdb_id_col, gdp_b_lst + gtp_1_lst),
        get_file_path(f"{cluster_col}_{pymol_pml_file}", dir_path=out_path),
        stick_resids=[32],
        loop_resids=[sw1_resids],
        group_col=nuc_class_col,
        color_palette=nuc_color_dict,
        thick_bb=False,
        show_bio=True,
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
            df,
            get_file_path(pymol_pml_file, dir_str=loop_name, dir_path=out_path),
            group_col=loop_name,
            stick_resids=stick_resids,
            loop_resids=loop_resids,
            style_ribbon=True,
            thick_bb=False,
            show_bio=True,
            color_palette=conf_color_dict[loop_name],
            sup_group=True,
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
                ),
                group_col=loop_name,
                stick_resids=stick_resids,
                loop_resids=loop_resids,
                style_ribbon=True,
                thick_bb=False,
                show_bio=True,
                show_prot=prot_color_dict[prot_class.split(".")[0]],
                color_palette=conf_color_dict[loop_name],
                sup_group=True,
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
                sup_group=True,
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
            sup_group=True,
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

    sw1_gtp_df = mask_unequal(df, sw1_name, [sw1_nf_name, sw1_gdp_name, noise_name])

    dist_table_path = get_file_path(
        dist_table_file, dir_str=y32_name, dir_path=out_path
    )

    build_dist_table(
        sw1_gtp_df,
        x_resids=[32],
        y_resids=[bio_lig_col],
        x_atomids=["OH"],
        y_atomids=[gtp_atomids],
        atom_dist_col_lst=[atom_dist_col],
        hb_status_col_lst=[hb_status_col],
        outlier_col_lst=[outlier_col],
        dist_table_path=dist_table_path,
        check_hb=True,
    )

    dist_df = load_table(dist_table_path)

    dist_df[hb_status_col] = dist_df[hb_status_col].map(sw1_gtp_dict)

    make_facet_plot(
        dist_df,
        get_file_path(plot_img_file, dir_str=y32_name, dir_path=out_path),
        x_col=atom_dist_col,
        x_str="Y32(OH):3P(O1G)\nDistance (Å)",
        plot_width=1.5,
        plot_height=1.5,
        plot_kde=True,
        show_legend=False,
        x_round=1,
        y_round=1,
        x_ticks=[2.5, 5.0, 7.5, 10.0],
        hue_col=hb_status_col,
        hue_palette=sw1_gtp_color_dict,
    )

    hras_df = mask_equal(
        mask_equal(mask_equal(dist_df, bio_lig_col, "GNP"), outlier_col, "False"),
        gene_class_col,
        hras_name,
    )
    write_pymol_script(
        hras_df,
        get_file_path(
            f"{hras_name}_{pymol_pml_file}", dir_str=y32_name, dir_path=out_path
        ),
        stick_resids=[32],
        group_col=hb_status_col,
        color_palette=sw1_gtp_color_dict,
        thick_bb=False,
        show_bio=True,
        sup_coord_path=sup_core_path,
        sup_chainid=sup_chainid,
        set_view=y32_view,
        x_hb_resids=[32],
        x_hb_atomids=["OH"],
        y_hb_resids=[bio_lig_col],
        y_hb_atomids=[gtp_atomids],
        show_hb=[sw1_gtp_dir_color],
        show_wmhb=[sw1_gtp_wat_color],
        show_resids=show_resids,
    )

    r_df = mask_equal(dist_df, sw2_name, sw2_gtp_r_name)

    dih_table_path = get_file_path(
        f"{q61_name}_{dih_table_file}",
        dir_str=y32_name,
        dir_path=out_path,
    )

    build_dih_table(
        r_df, dih_dict, chi1_resids=61, chi2_resids=61, dih_table_path=dih_table_path
    )

    dih_df = load_table(dih_table_path)

    for mut_class in ["WT", "G12D", "G12V"]:

        mut_df = mask_equal(
            dih_df,
            mut_class_col,
            mut_class,
        )

        make_facet_plot(
            mut_df,
            get_file_path(
                f"{mut_class}_{q61_name}_{chi1_col}_{chi2_col}_{plot_img_file}",
                dir_str=y32_name,
                dir_path=out_path,
            ),
            x_col=chi1_col,
            y_col=chi2_col,
            hue_col=hb_status_col,
            hue_palette=sw1_gtp_color_dict,
            x_str=f"χ$^1$",
            y_str=f"χ$^2$",
            plot_width=1.5,
            plot_height=1.5,
            kde_bw=1,
            show_legend=False,
        )

        write_pymol_script(
            mut_df,
            get_file_path(
                f"{mut_class}_{pymol_pml_file}", dir_str=y32_name, dir_path=out_path
            ),
            stick_resids=[12, 32, 61],
            group_col=hb_status_col,
            color_palette=sw1_gtp_color_dict,
            thick_bb=False,
            show_bio=True,
            sup_coord_path=sup_core_path,
            sup_chainid=sup_chainid,
            set_view=mut_view,
            x_hb_resids=[12, 12, 32, 61],
            x_hb_atomids=["OD1", "OD1", "OH", "OE1"],
            y_hb_resids=[32, bio_lig_col, bio_lig_col, 32],
            y_hb_atomids=["OH", gtp_atomids, gtp_atomids, "OH"],
            show_hb=[
                sw1_gtp_wat_color,
                sw1_gtp_wat_color,
                sw1_gtp_dir_color,
                sw1_gtp_dir_color,
            ],
            show_wmhb=[sw1_gtp_wat_color] * 4,
            show_resids=show_resids,
        )


def plot_pockets(df, pocket_df, dih_dict, sup_core_path, out_path):

    dih_table_path = get_file_path(
        f"{y71_name}_{dih_table_file}",
        dir_str=pocket_class_col,
        dir_path=out_path,
    )

    build_dih_table(pocket_df, dih_dict, chi1_resids=71, dih_table_path=dih_table_path)

    dih_df = load_table(dih_table_path)

    make_facet_plot(
        build_dih_table(
            df,
            dih_dict,
            chi1_resids=71,
        ),
        get_file_path(
            f"{nuc_class_col}_{y71_name}_{chi1_col}_{plot_img_file}",
            dir_path=out_path,
        ),
        x_col=chi1_col,
        x_str="Y71 χ$^1$",
        hue_col=nuc_class_col,
        hue_palette=nuc_color_dict,
        plot_width=1.5,
        plot_height=1.5,
        y_round=2,
        kde_bw=0.5,
        show_legend=True,
        legend_cols=1,
        legend_pad=3,
    )

    make_facet_plot(
        dih_df,
        get_file_path(
            f"{y71_name}_{chi1_col}_{plot_img_file}",
            dir_str=pocket_class_col,
            dir_path=out_path,
        ),
        x_col=chi1_col,
        x_str="Y71 χ$^1$",
        col_col=pocket_site_col,
        col_order=pocket_site_lst,
        col_wrap=3,
        hue_col=pocket_type_col,
        hue_palette=pocket_color_dict,
        hue_order=list(pocket_color_dict.keys()),
        plot_width=2.6,
        plot_height=1.2,
        y_round=2,
        kde_bw=0.5,
    )

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
                mask_equal(dih_df, pharm_class_col, pocket_site),
                get_file_path(
                    f"{pocket_site}_{y71_name}_{pymol_pml_file}",
                    dir_str=pocket_class_col,
                    dir_path=out_path,
                ),
                stick_resids=[71],
                loop_resids=[71],
                color_palette=[sw2_color],
                group_col=rotamer_col,
                color_group=False,
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
                sw2_gtp_r_name,
                sw2_gtp_sp12_a_name,
                sw2_gtp_sp12_b_name,
                sw2_gdp_sp12_name,
                noise_name,
            ]
        },
        sp2_name: {
            sw2_name: [
                sw2_gdp_sp2_a_name,
                sw2_gdp_sp2_b_name,
                noise_name,
            ]
        },
    }

    pocket_loop_cluster_dict[sp12_name][sw1_name] = [
        sw1_gdp_name,
        sw1_gtp_wat_name,
        noise_name,
    ]
    pocket_loop_cluster_dict[sp2_name][sw1_name] = [
        sw1_gdp_name,
        noise_name,
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

    # plot_pymol(df, sup_core_path, out_path)
    # plot_dist(df, out_path)
    plot_pockets(df, pocket_df, dih_dict, sup_core_path, out_path)
    # plot_rama(df, dih_dict, out_path,num_cpu=num_cpu)
