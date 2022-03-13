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

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from statannot import add_stat_annotation
import seaborn as sns

import warnings

warnings.filterwarnings("ignore", module="seaborn")

from ..functions.lst import type_lst, format_nested_lst
from ..functions.color import change_hex_alpha, get_lst_colors, gray_hex
from ..functions.col import (
    reformat_col_lst,
    resid_col,
    bb_col_lst,
    sc_col_lst,
    phi_col,
    psi_col,
)
from ..functions.table import (
    reformat_val_table,
    lst_col,
    mask_unequal,
    title_str,
    get_ncols,
)
from ..functions.plot import prep_plot_col
from ..functions.path import append_file_path

grid_hex = change_hex_alpha(gray_hex, 0.25)

bb_lim = (-180, 180)
bb_ticks = [-90, 0, 90]
sc_lim = (0, 360)
sc_ticks = [120, 240]

bb_grid_lst = [
    [[-180, 0], [50, 50], [180, 180]],
    [[-180, 0], [-180, -180], [-100, -100]],
    [[-180, 0], [-100, -100], [50, 50]],
    [[0, 180], [-180, -180], [-50, -50]],
    [[0, 180], [100, 100], [180, 180]],
    [[0, 180], [-50, -50], [100, 100]],
]

sc_x_lst = [[[0, 120]], [[120, 240]], [[240, 360]]]

sc_y_lst = [[[0, 0], [120, 120]], [[120, 120], [240, 240]], [[240, 240], [360, 360]]]

sc_grid_lst = list()
for sc_x in sc_x_lst:
    for sc_y in sc_y_lst:
        sc_grid_lst.append(sc_x + sc_y)


def make_facet_plot(
    plot_df,
    x_col,
    plot_path=None,
    y_col=None,
    rename_x=None,
    x_order=None,
    x_palette=None,
    x_count=True,
    hue_col=None,
    rename_hue=None,
    hue_order=None,
    hue_palette=None,
    hue_count=True,
    row_col=None,
    rename_row=None,
    row_order=None,
    row_count=True,
    row_palette=None,
    col_col=None,
    rename_col=None,
    col_order=None,
    col_count=False,
    col_wrap=None,
    col_palette=None,
    plot_width=4,
    plot_height=None,
    font_size=7,
    marker_size=1,
    line_width=0.5,
    tick_len=2,
    show_legend=False,
    legend_pad=5,
    legend_cols=None,
    legend_marker_size=None,
    color_legend_text=False,
    highlight_palette=None,
    highlight_lst=None,
    darken_lst=None,
    darken_palette=None,
    x_lim=None,
    y_lim=None,
    x_ticks=None,
    y_ticks=None,
    x_round=0,
    y_round=0,
    x_rotation=0,
    y_rotation=0,
    x_ha="center",
    y_ha="right",
    x_str=None,
    y_str=None,
    x_pad=None,
    y_pad=None,
    h_lines=None,
    v_lines=None,
    h_color=None,
    v_color=None,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
    plot_scatter=True,
    plot_line=False,
    plot_reg=False,
    plot_kde=False,
    plot_kind=None,
    log_reg=False,
    reg_estimator=None,
    reg_bins=None,
    trun_reg=False,
    kde_common_norm=True,
    kde_common_grid=False,
    kde_bw=1.0,
    kde_thresh=0.05,
    kde_levels=10,
    kde_cut=3.0,
    kde_alpha=0.1,
    rug_alpha=0.5,
    rug_height=0.05,
    rug_width=0.05,
    stat_pairs=None,
    stat_test="t-test_ind",
    stat_loc="inside",
    stat_format="star",
    tick_mult=0.75,
    all_ticks=False,
):

    df = plot_df.copy(deep=True)

    df_col_lst = list(df.columns)

    reformat_lst = list()
    if x_col in reformat_col_lst:
        if x_col not in df_col_lst:
            reformat_lst.append(x_col)
    if y_col in reformat_col_lst:
        if y_col is not None:
            if y_col not in df_col_lst:
                reformat_lst.append(y_col)

    if len(reformat_lst) > 0:
        df = reformat_val_table(df, reformat_lst)

    if plot_kind is None:
        df[x_col] = df[x_col].map(float)
    else:
        if type(x_palette) == str:
            x_palette = [x_palette] * len(lst_col(df, x_col))
        df, x_lst, x_palette = prep_plot_col(
            df,
            x_col,
            color_palette=x_palette,
            rename_vals=rename_x,
            order_lst=x_order,
            label_count=x_count,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
        )

    if y_col is not None:
        df[y_col] = df[y_col].map(float)

    if hue_col is not None:
        df, hue_lst, hue_palette = prep_plot_col(
            df,
            hue_col,
            color_palette=hue_palette,
            rename_vals=rename_hue,
            order_lst=hue_order,
            label_count=hue_count,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
        )
        if len(hue_lst) == 1:
            if type(show_legend) != dict and show_legend == False:
                show_legend = False
    else:
        hue_lst = type_lst(0)
        if type(show_legend) != dict:
            show_legend = False

    if row_col is not None:
        df, row_lst, row_palette = prep_plot_col(
            df,
            row_col,
            color_palette=row_palette,
            rename_vals=rename_row,
            order_lst=row_order,
            label_count=row_count,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
        )

        total_row = len(row_lst)

        if row_palette is not None:
            row_color_dict = get_lst_colors(
                row_lst, palette=row_palette, return_dict=True
            )

        if hue_col == row_col:
            if type(show_legend) != dict:
                show_legend = False
    else:
        row_lst = None
        total_row = 1

    title_col = True
    if col_col is None:
        col_col = resid_col
        if resid_col not in list(df.columns):
            df[resid_col] = "None"
            title_col = False

    df, col_lst, col_palette = prep_plot_col(
        df,
        col_col,
        color_palette=col_palette,
        rename_vals=rename_col,
        order_lst=col_order,
        label_count=col_count,
        count_chain=count_chain,
        count_pdb=count_pdb,
        count_cf=count_cf,
    )
    total_col = len(col_lst)

    if col_palette is not None:
        col_color_dict = get_lst_colors(col_lst, palette=col_palette, return_dict=True)

    df = mask_unequal(df, x_col, 999.00)
    if x_str is None:
        x_str = title_str(x_col)
    if x_col in bb_col_lst:
        x_lim = bb_lim
        if x_ticks is not False:
            x_ticks = bb_ticks
    elif x_col in sc_col_lst and x_lim is None:
        x_lim = sc_lim
        if x_ticks is not False:
            x_ticks = sc_ticks
    if y_col is not None:
        df = mask_unequal(df, y_col, 999.00)
        if y_str is None:
            y_str = title_str(y_col)
        if y_col in bb_col_lst:
            y_lim = bb_lim
            if y_ticks is not False:
                y_ticks = bb_ticks
        elif y_col in sc_col_lst and y_lim is None:
            y_lim = sc_lim
            if y_ticks is not False:
                y_ticks = sc_ticks
    else:
        if y_str is None:
            y_str = "Density"

    if x_lim is not None:
        if type(x_lim) != tuple:
            x_lim = tuple(x_lim)
    if y_lim is not None:
        if type(y_lim) != tuple:
            y_lim = tuple(y_lim)

    if plot_kind is None or hue_col is not None:
        hue_color_dict = get_lst_colors(hue_lst, palette=hue_palette, return_dict=True)
    else:
        hue_color_dict = get_lst_colors(x_lst, palette=x_palette, return_dict=True)
    if highlight_lst is not None:
        highlight_lst = format_nested_lst(highlight_lst, return_str=True)
        highlight_color_lst = get_lst_colors(highlight_lst, palette=highlight_palette)

    if darken_lst is not None:
        darken_lst = format_nested_lst(darken_lst, return_str=True)
        darken_color_lst = get_lst_colors(darken_lst, palette=darken_palette)

    sns.set_context("paper")
    sns.set_style("ticks")
    sns.set_palette(list(hue_color_dict.values()))

    bbox_extra_artists = tuple()

    if row_col is None:
        if col_wrap is None:
            col_wrap = total_col
        total_col = col_wrap
        total_row = total_col / col_wrap
        show_margins = False
    else:
        col_wrap = None
        show_margins = True

    if plot_kind is None:

        g = sns.FacetGrid(
            df,
            hue=hue_col,
            hue_order=hue_lst,
            row=row_col,
            row_order=row_lst,
            col=col_col,
            col_order=col_lst,
            col_wrap=col_wrap,
            margin_titles=show_margins,
        )

        if y_col is None:
            g.map(
                sns.distplot,
                x_col,
                kde=True,
                hist=False,
                rug=True,
                kde_kws={
                    "common_norm": kde_common_norm,
                    "common_grid": kde_common_grid,
                    "bw_adjust": kde_bw,
                    "thresh": kde_thresh,
                    "levels": kde_levels,
                    "cut": kde_cut,
                    "alpha": kde_alpha,
                    "fill": True,
                    "linewidth": line_width,
                },
                rug_kws={"lw": rug_width, "alpha": rug_alpha, "height": rug_height},
            )
        else:

            if plot_line:
                g.map(sns.lineplot, x_col, y_col)
            if plot_reg:
                g.map(
                    sns.regplot,
                    x_col,
                    y_col,
                    scatter_kws={
                        "s": marker_size,
                        "linewidth": 0,
                        "alpha": 0.75,
                    },
                    line_kws={"linewidth": line_width},
                    logx=log_reg,
                    x_estimator=reg_estimator,
                    x_bins=reg_bins,
                    truncate=trun_reg,
                )
            if plot_kde:
                g.map(
                    sns.kdeplot,
                    x_col,
                    y_col,
                    common_norm=kde_common_norm,
                    common_grid=kde_common_grid,
                    bw_adjust=kde_bw,
                    cut=kde_cut,
                    thresh=kde_thresh,
                    levels=kde_levels,
                    fill=True,
                    linewidth=0,
                    alpha=0.25,
                )

            if plot_scatter:
                g.map(
                    plt.scatter,
                    x_col,
                    y_col,
                    s=marker_size,
                    linewidth=0,
                    alpha=0.75,
                )
    else:
        if plot_kind == "strip" or plot_kind == "swarm" or plot_kind == "point":
            g = sns.catplot(
                data=df,
                x=x_col,
                y=y_col,
                hue=hue_col,
                row=row_col,
                col=col_col,
                col_wrap=col_wrap,
                order=x_lst,
                hue_order=hue_lst,
                row_order=row_lst,
                col_order=col_lst,
                kind=plot_kind,
                palette=hue_color_dict,
                legend=False,
                sharex=False,
                sharey=True,
                margin_titles=True,
                s=marker_size,
                linewidth=0,
                alpha=0.75,
                dodge=True,
            )
        else:
            if (
                row_col is not None
                or col_col is not resid_col
                or plot_kind == "count"
                or plot_kind == "bar"
            ):

                if plot_kind == "count":
                    cat_width = 0
                else:
                    cat_width = line_width

                g = sns.catplot(
                    data=df,
                    x=x_col,
                    y=y_col,
                    hue=hue_col,
                    row=row_col,
                    col=col_col,
                    col_wrap=col_wrap,
                    order=x_lst,
                    hue_order=hue_lst,
                    row_order=row_lst,
                    col_order=col_lst,
                    kind=plot_kind,
                    palette=hue_color_dict,
                    legend=False,
                    sharex=False,
                    sharey=True,
                    linewidth=cat_width,
                    margin_titles=True,
                )
            else:
                if plot_kind == "box" or plot_kind == "boxen" or plot_kind == "violin":

                    g = sns.catplot(
                        data=df,
                        x=x_col,
                        y=y_col,
                        hue=hue_col,
                        order=x_lst,
                        hue_order=hue_lst,
                        kind=plot_kind,
                        color="white",
                        legend=False,
                        linewidth=line_width,
                        showfliers=False,
                    )

                    ax = sns.stripplot(
                        data=df,
                        x=x_col,
                        y=y_col,
                        hue=hue_col,
                        order=x_lst,
                        hue_order=hue_lst,
                        palette=hue_color_dict,
                        s=marker_size,
                        alpha=0.75,
                        split=True,
                    )
                    if ax.get_legend() is not None:
                        ax.get_legend().remove()

                    if stat_pairs is not None:
                        x_dict = dict()
                        for val in x_lst:
                            x_dict[val.split(" (N=")[0]] = val

                        box_pairs = list()
                        for i, stat in enumerate(stat_pairs):
                            box_pairs.append((x_dict[stat[0]], x_dict[stat[1]]))

                        add_stat_annotation(
                            ax,
                            data=df,
                            x=x_col,
                            y=y_col,
                            order=x_lst,
                            box_pairs=box_pairs,
                            test=stat_test,
                            text_format=stat_format,
                            loc=stat_loc,
                            verbose=2,
                            fontsize=font_size,
                            linewidth=line_width,
                        )

    if plot_height is None:
        plot_height = plot_width * (total_row / total_col)

    g.fig.set_figheight(plot_height)
    g.fig.set_figwidth(plot_width)

    if plot_kind is None:
        if (x_ticks is None or x_ticks is False) and (
            y_ticks is None or y_ticks is False
        ):
            g.set(
                xlim=x_lim,
                ylim=y_lim,
            )
        elif x_ticks is None or x_ticks is False:
            g.set(
                xlim=x_lim,
                ylim=y_lim,
                yticks=y_ticks,
            )
        elif y_ticks is None or y_ticks is False:
            g.set(
                xlim=x_lim,
                ylim=y_lim,
                xticks=x_ticks,
            )
        else:
            g.set(
                xlim=x_lim,
                ylim=y_lim,
                xticks=x_ticks,
                yticks=y_ticks,
            )
    else:
        if y_ticks is not None and y_lim is not None:
            g.set(
                ylim=y_lim,
                yticks=y_ticks,
            )
        elif y_ticks is not None:
            g.set(
                yticks=y_ticks,
            )
        elif y_lim is not None:
            g.set(
                ylim=y_lim,
            )

    g.fig.subplots_adjust(wspace=0, hspace=0)

    g.fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)

    x_label = plt.xlabel(x_str, fontsize=font_size, labelpad=x_pad)
    y_label = plt.ylabel(y_str, fontsize=font_size, labelpad=y_pad)

    bbox_extra_artists += (x_label, y_label)

    if show_legend or type(show_legend) == dict:

        if legend_cols is None:
            legend_cols = get_ncols(hue_lst)

        if legend_marker_size is None:
            legend_marker_size = marker_size * 3

        if type(show_legend) == dict:
            handles = [
                plt.plot(
                    [],
                    [],
                    marker="o",
                    ls="",
                    markersize=legend_marker_size,
                    markeredgewidth=0,
                    markerfacecolor=show_legend[hue],
                    label=hue,
                    color=show_legend[hue],
                )[0]
                for hue in list(show_legend.keys())
            ]
        else:
            handles = [
                plt.plot(
                    [],
                    [],
                    marker="o",
                    ls="",
                    markersize=legend_marker_size,
                    markeredgewidth=0,
                    markerfacecolor=hue_color_dict[hue],
                    label=hue,
                    color=hue_color_dict[hue],
                )[0]
                for hue in hue_lst
            ]

        legend_label_color = None
        if color_legend_text:
            legend_label_color = "linecolor"

        legend = g.fig.legend(
            handles=handles,
            fontsize=font_size,
            ncol=legend_cols,
            loc="upper center",
            frameon=False,
            bbox_transform=g.fig.transFigure,
            bbox_to_anchor=(0.5, 0),
            borderaxespad=legend_pad,
            labelcolor=legend_label_color,
        )

        bbox_extra_artists += (legend,)

    g.set_axis_labels("", "")

    tick_index = 0

    total_axes = len(col_lst)
    if row_col is not None:
        total_axes *= len(row_lst)

    ax_lst = col_lst
    if row_col is not None and col_col is None:
        ax_lst = row_lst
    elif row_col is not None and col_col is not None:
        ax_lst = list()
        for row in row_lst:
            for col in col_lst:
                ax_lst.append((row, col))

    if h_lines is not None:
        h_line_lst = type_lst(h_lines)
        if h_color is None:
            h_color = "black"

    if v_lines is not None:
        v_line_lst = type_lst(v_lines)
        if v_color is None:
            v_color = "black"

    for i, ax in enumerate(g.axes.flat):

        ax_name = ax_lst[i]

        if h_lines is not None:
            for h in h_line_lst:
                ax.axhline(y=h, linewidth=line_width, linestyle="--", color=h_color)

        if v_lines is not None:
            for v in v_line_lst:
                ax.axvline(x=v, linewidth=line_width, linestyle="--", color=v_color)

        clean_name = ax_name
        if type(clean_name) == str:
            if " (N=" in clean_name:
                clean_name = clean_name.split(" (N=")[0]
        elif type(clean_name) == tuple:
            clean_name = list(clean_name)
            for n, name in enumerate(clean_name):
                clean_name[n] = name.split(" (N=")[0]
            clean_name = tuple(clean_name)

        ax_grid_hex = grid_hex
        if darken_lst is not None:
            for d, darken in enumerate(darken_lst):
                if clean_name in darken:
                    ax.set_facecolor(darken_color_lst[d])
                    ax_grid_hex = "white"

        if x_col == phi_col and y_col == psi_col:
            for bb_grid in bb_grid_lst:
                ax.fill_between(
                    bb_grid[0],
                    bb_grid[1],
                    bb_grid[2],
                    facecolor="none",
                    edgecolor=ax_grid_hex,
                    linewidth=line_width,
                    zorder=0,
                )
        elif (x_col in sc_col_lst and y_col in sc_col_lst) and (
            x_ticks is None and y_ticks is None
        ):
            for sc_grid in sc_grid_lst:
                ax.fill_between(
                    sc_grid[0],
                    sc_grid[1],
                    sc_grid[2],
                    facecolor="none",
                    edgecolor=ax_grid_hex,
                    linewidth=line_width,
                    zorder=0,
                )
        else:
            if plot_kind is None:
                ax.xaxis.grid(color=ax_grid_hex, linewidth=line_width)
            ax.yaxis.grid(color=ax_grid_hex, linewidth=line_width)

        ax.set_axisbelow(True)

        if row_col is None:
            if ax.get_title():
                ax.set_title("")

            if len(col_lst) > 1:

                highlight_color = gray_hex
                if highlight_lst is not None:
                    for h, highlight in enumerate(highlight_lst):
                        if clean_name in highlight:
                            highlight_color = highlight_color_lst[h]

                divider = make_axes_locatable(ax)
                cax = divider.append_axes("top", size="20%", pad=0)
                cax.get_xaxis().set_visible(False)
                cax.get_yaxis().set_visible(False)
                cax.set_facecolor(highlight_color)
                for _, spine in cax.spines.items():
                    spine.set_visible(True)
                    spine.set_color("black")
                    spine.set_linewidth(line_width)
                cax.text(
                    0.5,
                    0.4,
                    str(ax_name),
                    color="white",
                    fontsize=font_size * 0.75,
                    verticalalignment="center",
                    horizontalalignment="center",
                    transform=cax.transAxes,
                )
        else:
            if ax.get_title():

                if title_col:
                    col = ax.get_title().split("=", 1)[1]
                else:
                    col = ""

                col_color = "black"
                if col_palette is not None:
                    col_color = col_color_dict[col]

                col_txt = ax.set_title(
                    col,
                    fontsize=font_size,
                    color=col_color,
                )

                bbox_extra_artists += (col_txt,)

            if ax.texts:
                row = ax.texts[0].get_text().split("=", 1)[1]

                row_color = "black"
                if row_palette is not None:
                    row_color = row_color_dict[row[1:]]

                row_txt = ax.text(
                    ax.texts[0].get_unitless_position()[0],
                    ax.texts[0].get_unitless_position()[1],
                    row,
                    transform=ax.transAxes,
                    va="center",
                    fontsize=font_size,
                    color=row_color,
                )
                ax.texts[0].remove()

                bbox_extra_artists += (row_txt,)

        add_ticks = False
        bottom_label = True
        left_label = True
        bottom = True
        left = True

        if total_col == 1:
            add_ticks = True
            if i != total_axes - 1:
                bottom_label = False
                bottom = False
        elif total_axes == 1 and row_col is None:
            add_ticks = True
        elif total_axes == col_wrap and row_col is None:
            add_ticks = True
            if i != 0:
                left_label = False
                left = False
        elif i % total_col == 0:
            if tick_index % 2 != 0 or all_ticks:
                add_ticks = True
            tick_index += 1
            if not all_ticks:
                bottom_label = False
                bottom = False
        elif (total_axes - i) < total_col + 1:
            if i % 2 != 0 or all_ticks:
                add_ticks = True
            left_label = False
            left = False
        else:
            bottom_label = False
            left_label = False
            bottom = False
            left = False

        if add_ticks:

            tick_format = "{:,."

            if plot_kind is None:
                x_tick_format = tick_format + str(x_round)

                x_tick_format += "f}"

                x_tick_lst = ax.get_xticks()
                ax.xaxis.set_major_locator(mticker.FixedLocator(x_tick_lst))
                ax.set_xticklabels(
                    [x_tick_format.format(x) for x in x_tick_lst],
                    fontsize=font_size * tick_mult,
                    rotation=x_rotation,
                    ha=x_ha,
                )
            else:
                ax.set_xticklabels(
                    ax.get_xticklabels(),
                    fontsize=font_size * tick_mult,
                    rotation=45,
                    ha="right",
                )

            y_tick_format = tick_format + str(y_round)
            y_tick_format += "f}"

            y_tick_lst = ax.get_yticks()
            ax.yaxis.set_major_locator(mticker.FixedLocator(y_tick_lst))
            ax.set_yticklabels(
                [y_tick_format.format(y) for y in y_tick_lst],
                fontsize=font_size * tick_mult,
                rotation=y_rotation,
                ha=y_ha,
            )

        elif not add_ticks:
            bottom_label = False
            left_label = False

        if x_ticks is False:
            bottom_label = False
            bottom = False
        if y_ticks is False:
            left = False
            left_label = False

        ax.tick_params(
            direction="out",
            labelbottom=bottom_label,
            labelleft=left_label,
            bottom=bottom,
            left=left,
            length=tick_len,
            width=line_width,
            colors="black",
        )

        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(line_width)

    if plot_path is None:
        return g
    else:
        append_file_path(plot_path)

        if 'png' in plot_path:
            plot_format = 'png'
        if 'pdf' in plot_path:
            plot_format = 'pdf'

        plt.savefig(
            plot_path,
            format=plot_format,
            bbox_extra_artists=bbox_extra_artists,
            bbox_inches="tight",
            pad_inches=0.0,
            dpi=600,
        )
        plt.close()

        print("Made facet plot!")
