# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the pdbcore project.

The pdbcore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import seaborn as sns

from util.plot import prep_plot_col
from util.path import save_table
from util.data import (
    build_col_count_dict,
    get_ncols,
    convert_col_percent,
    mask_greater,
    mask_unequal,
)
from util.path import append_file_path
from util.lst import type_lst
from util.sig import calc_rr
from util.col import pdb_id_col, total_col, index_col, risk_ratio_col, sig_col


def make_heatmap_plot(
    df,
    plot_path,
    x_col=None,
    y_col=None,
    x_order=None,
    rename_x=None,
    x_count=False,
    y_order=None,
    rename_y=None,
    y_count=False,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
    plot_width=4,
    plot_height=4,
    font_size=10,
    line_width=0.5,
    tick_len=2,
    cmap_palette="viridis",
    x_str=None,
    y_str=None,
    v_min=None,
    v_max=None,
    x_ticks="auto",
    y_ticks="auto",
    col_percent=False,
    row_percent=False,
    add_margins=False,
    show_annot=False,
    annot_round=None,
    row_colors=None,
    col_colors=None,
    row_labels=None,
    col_labels=None,
    colors_ratio=0.05,
    show_cbar=True,
    cbar_x=0.2,
    cbar_y=0,
    cbar_len=0.4,
    cbar_thick=0.03,
    cbar_extend="both",
    cbar_str=None,
    legend_dict=None,
    legend_cols=None,
    legend_pad=1,
    legend_x=0.5,
    legend_y=1,
    show_sig=False,
    min_rr=1,
    exp_col=None,
    out_col=None,
    h_line_lst=None,
    v_line_lst=None,
    sig_table_path=None,
    show_grid=False,
    grid_color="white",
    add_border=False,
    highlight_lst=None,
    highlight_color="red",
    id_col=None,
):
    if x_col is not None and y_col is not None:

        if id_col is None:
            df = df.reset_index()
            id_col = index_col

        df, x_order = prep_plot_col(
            df,
            x_col,
            rename_vals=rename_x,
            order_lst=x_order,
            label_count=x_count,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
            return_palette=False,
        )
        df, y_order = prep_plot_col(
            df,
            y_col,
            rename_vals=rename_y,
            order_lst=y_order,
            label_count=y_count,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
            return_palette=False,
        )

        if col_percent or row_percent:
            x_count_dict = build_col_count_dict(df, x_col, col_lst=id_col)
            y_count_dict = build_col_count_dict(df, y_col, col_lst=id_col)

        df = pd.pivot_table(
            df, index=y_col, columns=x_col, values=id_col, aggfunc="nunique"
        ).fillna(0)

        if col_percent or row_percent:
            for index in list(df.index.values):
                for col in list(df.columns):
                    if col_percent:
                        total = x_count_dict[col]
                    elif row_percent:
                        total = y_count_dict[index]
                    df.at[index, col] = round((df.at[index, col] / total) * 100, 1)

        for x in x_order:
            if x not in list(df.columns):
                for col in list(df.columns):
                    df[x] = np.nan

        for y in y_order:
            if y not in list(df.index.values):
                for col in list(df.columns):
                    df.at[y, col] = np.nan

        if x_order is not None and y_order is not None:
            df = df.loc[y_order, x_order]
        elif x_order is not None and y_order is None:
            df = df.loc[:, x_order]
        elif x_order is None and y_order is not None:
            df = df.loc[y_order, :]

    for col in list(df.columns):
        df[col] = df[col].map(float)

    if show_sig:

        if x_col is None:
            x_col = str(0)
        if y_col is None:
            y_col = str(1)

        if exp_col is None:
            exp_col = x_col
        if out_col is None:
            out_col = y_col

        stat_df = df.copy(deep=True)

        stat_df = stat_df.reset_index()

        stat_df = stat_df.rename(columns={index_col: y_col})

        stat_df = pd.melt(
            stat_df,
            id_vars=y_col,
            var_name=x_col,
            value_name=total_col,
        )

        stat_df = calc_rr(stat_df, exp_col, out_col, correct_method="fdr_bh")
        stat_df = mask_greater(stat_df, risk_ratio_col, min_rr)
        stat_df = mask_unequal(stat_df, sig_col, "ns")

        df_index_lst = list(df.index.values)
        df_col_lst = list(df.columns)

        stat_df = stat_df.reset_index(drop=True)

        highlight_lst = list()
        for index in list(stat_df.index.values):
            highlight_lst.append(
                (
                    df_index_lst.index(stat_df.at[index, y_col]),
                    df_col_lst.index(stat_df.at[index, x_col]),
                )
            )

        if sig_table_path is not None:
            save_table(sig_table_path, stat_df)

    if row_labels is not None:
        row_dict = dict()
        for color, label in zip(row_colors, row_labels):
            row_dict[label] = color
        row_colors = pd.DataFrame(row_dict, index=df.index.values)

    if col_labels is not None:
        col_dict = dict()
        for color, label in zip(col_colors, col_labels):
            col_dict[label] = color
        col_colors = pd.DataFrame(col_dict, index=df.columns)

    sns.set_context("paper")
    sns.set_style("ticks")

    if annot_round is None:
        fmt = "g"
    elif type(annot_round) == str:
        fmt = annot_round
    else:
        fmt = f"0.{annot_round}f"
        if type(show_annot) == bool:
            if show_annot:
                show_annot = df.copy(deep=True)
                for index in list(df.index.values):
                    for col in list(df.columns):
                        show_annot.at[index, col] = str(
                            round(show_annot.at[index, col], annot_round)
                        ).replace("." + "0" * annot_round, "")
                show_annot = show_annot.fillna("-")

    if type(show_annot) != bool:
        show_annot = show_annot
        fmt = ""

    if x_col is None and y_col is None:
        if col_percent:
            for col in list(df.columns):
                df = convert_col_percent(df, col)
                df[col] = df[col].round(1)
            if add_margins:
                df.loc["All"] = 100
            if row_colors is not None:
                row_colors.append("white")
        elif row_percent:
            df = df.div(df.sum(axis=1), axis=0) * 100
            for col in list(df.columns):
                df[col] = df[col].round(1)
            if add_margins:
                df["All"] = 100
            if col_colors is not None:
                col_colors.append("white")

    grid_width = 0
    if show_grid:
        grid_width = line_width

    df.replace({int(0): np.nan}, inplace=True)
    df.replace({float(0): np.nan}, inplace=True)

    g = sns.clustermap(
        df,
        figsize=(plot_width, plot_height),
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        col_colors=col_colors,
        colors_ratio=colors_ratio,
        cbar_pos=(cbar_x, cbar_y, cbar_len, cbar_thick),
        cmap=cmap_palette,
        annot_kws={"fontsize": font_size * 0.75},
        cbar_kws={
            "orientation": "horizontal",
            "extend": cbar_extend,
            "spacing": "uniform",
        },
        vmin=v_min,
        vmax=v_max,
        xticklabels=x_ticks,
        yticklabels=y_ticks,
        annot=show_annot,
        fmt=fmt,
        linewidths=grid_width,
        linecolor=grid_color,
    )

    if add_border:
        for _, spine in g.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(line_width)

    if highlight_lst is not None:
        highlight_lst = type_lst(highlight_lst)

        max_zorder = max([_.zorder for _ in g.ax_heatmap.get_children()])

        for highlight in highlight_lst:
            g.ax_heatmap.add_patch(
                Rectangle(
                    (highlight[1], highlight[0]),
                    1,
                    1,
                    edgecolor=highlight_color,
                    fill=False,
                    lw=line_width * 2,
                    zorder=max_zorder,
                )
            )

    if x_str is None:
        x_str = ""

    x_label = g.ax_heatmap.set_xlabel(x_str, fontsize=font_size)

    if y_str is None:
        y_str = ""

    y_label = g.ax_heatmap.set_ylabel(
        y_str, fontsize=font_size, rotation=-90, labelpad=10
    )

    bbox_extra_artists = (
        x_label,
        y_label,
    )

    if h_line_lst is not None:
        h_line_lst = type_lst(h_line_lst)
        index_lst = list(df.index.values)
        for h_index in [index_lst.index(x) for x in h_line_lst]:
            g.ax_heatmap.axhline(
                h_index,
                linewidth=line_width,
                color="black",
            )
    if v_line_lst is not None:
        v_line_lst = type_lst(v_line_lst)
        col_lst = list(df.columns)
        for v_index in [col_lst.index(x) for x in v_line_lst]:
            g.ax_heatmap.axvline(v_index, linewidth=line_width, color="black")

    plt.setp(
        g.ax_heatmap.get_xticklabels(),
        fontsize=font_size * 0.75,
        rotation=45,
        ha="right",
    )
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=font_size * 0.75)

    plt.tick_params(direction="out", length=tick_len, width=line_width, colors="black")

    if row_labels is not None:
        plt.setp(
            g.ax_row_colors.get_xticklabels(),
            fontsize=font_size * 0.75,
            rotation=45,
            ha="right",
        )

    if col_labels is not None:
        plt.setp(
            g.ax_col_colors.get_yticklabels(),
            fontsize=font_size * 0.75,
        )

    if show_cbar:

        cax = plt.gcf().axes[-1]
        cax.tick_params(
            labelsize=font_size * 0.75,
            length=tick_len,
            width=line_width,
            colors="black",
        )

        if cbar_str is not None:
            cax.set_xlabel(cbar_str, fontsize=font_size)

    else:
        g.cax.set_visible(False)

    if show_sig:
        sig_patch = mpatches.Patch(
            facecolor="none",
            edgecolor="red",
            linewidth=line_width,
            label=f"RR > {min_rr}; P < 0.05",
        )
        legend_1 = g.ax_col_dendrogram.legend(
            handles=[sig_patch],
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(legend_x, legend_y),
            fontsize=font_size * 0.75,
            borderaxespad=-legend_pad,
        )
        bbox_extra_artists += (legend_1,)

    elif legend_dict is not None:

        if legend_cols is None:
            legend_cols = get_ncols(list(legend_dict.keys()))

        for label, color in legend_dict.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
            legend_2 = g.ax_col_dendrogram.legend(
                ncol=legend_cols,
                frameon=False,
                loc="upper center",
                labelcolor="linecolor",
                bbox_to_anchor=(legend_x, legend_y),
                fontsize=font_size * 0.75,
                borderaxespad=-legend_pad,
            )
            bbox_extra_artists += (legend_2,)
    else:
        g.ax_col_dendrogram.set_visible(False)

    g.ax_row_dendrogram.set_visible(False)

    append_file_path(plot_path)

    plt.savefig(
        plot_path,
        format="pdf",
        bbox_extra_artists=bbox_extra_artists,
        bbox_inches="tight",
        pad_inches=0,
        dpi=600,
    )
    plt.close()

    print("Made heatmap plot!")
