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
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from .color import change_hex_alpha, get_lst_colors, gray_hex, black_hex
from .path import append_file_path
from .table import (
    format_val,
    make_dict,
    lst_col,
    build_label_dict,
    mask_equal,
    get_ncols,
)
from .col import pdb_id_col

grid_hex = change_hex_alpha(gray_hex, 0.25)


def prep_plot_col(
    df,
    col,
    color_palette=None,
    rename_vals=None,
    order_lst=None,
    label_count=False,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
    return_palette=True,
):

    df[col] = df[col].map(str)

    for index in list(df.index.values):

        df.at[index, col] = str(df.at[index, col]).split(" (N=")[0]

    if rename_vals is not None:
        if type(rename_vals) == dict:
            rename_dict = rename_vals
        elif type(rename_vals) == list:
            val_lst = lst_col(df, col, unique=True, return_str=True)
            rename_lst = format_val(rename_vals, return_str=True)
            rename_dict = make_dict(val_lst, rename_lst)
        df[col] = df[col].map(rename_dict)

    if order_lst is not None:
        order_lst = format_val(order_lst, return_str=True)
        df = mask_equal(df, col, order_lst)
        return_lst = order_lst.copy()
    else:
        return_lst = lst_col(df, col, unique=True, return_str=True)

    if label_count:
        col_dict = build_label_dict(
            df,
            col,
            return_str=True,
            count_chain=count_chain,
            count_pdb=count_pdb,
            count_cf=count_cf,
        )

        df[col] = df[col].map(col_dict)
        for i, row in enumerate(return_lst):
            return_lst[i] = col_dict[row]

    if type(color_palette) == dict:
        return_dict = dict()
        for val in return_lst:
            return_dict[val] = color_palette[val.split(" (")[0]]
    else:
        return_dict = color_palette

    if return_palette:
        return df, return_lst, return_dict
    else:
        return (
            df,
            return_lst,
        )


def make_legend_plot(
    plot_path,
    legend_dict,
    plot_width=2,
    plot_height=2,
    font_size=7,
    marker_shape="s",
    marker_size=3,
    legend_cols=None,
    legend_title=None,
    color_text=False,
):

    fig, ax = plt.subplots()

    fig.set_size_inches(plot_width, plot_height)

    handles = [
        plt.plot(
            [],
            [],
            marker=marker_shape,
            ls="",
            markersize=marker_size,
            markeredgewidth=0,
            markerfacecolor=color,
            label=label,
            color=color,
        )[0]
        for label, color in legend_dict.items()
    ]

    if legend_cols is None:
        legend_cols = get_ncols(list(legend_dict.keys()))

    label_color = None
    if color_text:
        label_color = "linecolor"

    legend = ax.legend(
        handles=handles,
        fontsize=font_size,
        ncol=legend_cols,
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.5),
        title=legend_title,
        labelcolor=label_color,
    )

    if legend_title is not None:
        plt.setp(legend.get_title(), fontsize=font_size)

    bbox_extra_artists = (legend,)

    plt.axis("off")

    append_file_path(plot_path)

    plt.savefig(
        plot_path,
        format="pdf",
        bbox_extra_artists=bbox_extra_artists,
        bbox_inches="tight",
        pad_inches=0.0,
        dpi=600,
    )
    plt.close()


def make_venn_plot(
    lst_1,
    lst_2,
    plot_path=None,
    label_1=None,
    label_2=None,
    color_1=None,
    color_2=None,
    color_inter=None,
    count_color="black",
    plot_title=None,
    plot_height=2,
    plot_width=2,
    font_size=7,
    alpha=0.75,
):

    if color_1 is None:
        color_1 = gray_hex
    if color_2 is None:
        color_2 = gray_hex
    if color_inter is None:
        color_inter = black_hex

    fig, ax = plt.subplots()

    fig.set_size_inches(plot_width, plot_height)

    total = len(set(lst_1 + lst_2))

    v = venn2(
        [set(lst_1), set(lst_2)],
        set_labels=(label_1, label_2),
        set_colors=(color_1, color_2),
        alpha=alpha,
        subset_label_formatter=lambda x: f"{x}\n({(x/total):1.0%})",
        ax=ax,
    )

    bbox_extra_artists = tuple()

    for text in v.set_labels:
        text.set_fontsize(font_size)
        bbox_extra_artists += (text,)

    for text in v.subset_labels:

        if text is not None:
            if text.get_text() == "0\n(0%)":
                text.set_text("")
            else:
                text.set_fontsize(font_size * 0.75)
                text.set_color(count_color)
                bbox_extra_artists += (text,)

    if plot_title is not None:
        title = fig.suptitle(plot_title, fontsize=font_size)
        bbox_extra_artists += (title,)

    if plot_path is None:
        return fig
    else:
        append_file_path(plot_path)

        plt.savefig(
            plot_path,
            format="pdf",
            bbox_extra_artists=bbox_extra_artists,
            bbox_inches="tight",
            pad_inches=0.0,
            dpi=600,
        )
        plt.close()


def make_stacked_barplot(
    plot_df,
    col_col,
    hue_col,
    plot_path,
    col_order=None,
    rename_col=None,
    col_count=False,
    hue_order=None,
    rename_hue=None,
    hue_count=False,
    hue_palette=None,
    x_str=None,
    y_str=None,
    font_size=7,
    plot_height=2,
    plot_width=2,
    line_width=0.5,
    show_legend=True,
    legend_pad=10,
    legend_cols=None,
    bar_width=0.5,
    bar_alpha=1,
    count_chain=True,
    count_pdb=False,
    count_cf=False,
    id_col=None,
    show_barh=False,
):

    df = plot_df.copy(deep=True)

    if id_col is None:
        id_column = pdb_id_col

    df, col_lst = prep_plot_col(
        df,
        col_col,
        rename_vals=rename_col,
        order_lst=col_order,
        label_count=col_count,
        count_chain=count_chain,
        count_pdb=count_pdb,
        count_cf=count_cf,
        return_palette=False,
    )

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

    hue_color_dict = get_lst_colors(hue_lst, palette=hue_palette, return_dict=True)

    sns.set_context("paper")
    sns.set_style("whitegrid")
    sns.set_palette(list(hue_color_dict.values()))

    df = pd.pivot_table(
        df,
        index=col_col,
        columns=hue_col,
        values=id_column,
        aggfunc="nunique",
    ).fillna(0)

    df = df.reindex(index=col_lst)
    df = df.reindex(columns=hue_lst)

    if show_barh:
        plot_kind = "barh"
        grid_axis = "y"
        if x_str is None:
            x_str = "% Structures"
        if y_str is None:
            y_str = col_col
    else:
        plot_kind = "bar"
        grid_axis = "x"
        if x_str is None:
            x_str = col_col
        if y_str is None:
            y_str = "% Structures"

    df.apply(lambda x: x / sum(x) * 100, axis=1).plot(
        kind=plot_kind,
        stacked=True,
        linewidth=0,
        width=bar_width,
        alpha=bar_alpha,
        figsize=(plot_width, plot_height),
        legend=show_legend,
    )

    if show_barh:
        plt.xticks(fontsize=font_size * 0.75)
    else:
        plt.xticks(fontsize=font_size * 0.75, rotation=45, ha="right")

    plt.yticks(fontsize=font_size * 0.75)
    plt.grid(axis=grid_axis, color=grid_hex, linewidth=line_width)

    x_label = plt.xlabel(x_str, fontsize=font_size)
    y_label = plt.ylabel(y_str, fontsize=font_size)

    bbox_extra_artists = (x_label, y_label)

    if show_legend or type(show_legend) == dict:

        if legend_cols is None:
            legend_cols = get_ncols(hue_lst)

        if type(show_legend) == dict:
            handles = [
                plt.plot(
                    [],
                    [],
                    marker="o",
                    ls="",
                    markeredgewidth=0,
                    markersize=3,
                    markerfacecolor=show_legend[hue],
                    label=hue,
                    color=hue,
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
                    markeredgewidth=0,
                    markersize=3,
                    markerfacecolor=hue_color_dict[hue],
                    label=hue,
                )[0]
                for hue in hue_lst
            ]

        legend = plt.legend(
            handles=handles,
            fontsize=font_size * 0.75,
            ncol=legend_cols,
            loc="upper center",
            frameon=False,
            bbox_to_anchor=(0.5, 0),
            borderaxespad=legend_pad,
        )

        bbox_extra_artists += (legend,)

    sns.despine(left=True)

    plt.savefig(
        plot_path,
        format="pdf",
        bbox_inches="tight",
        bbox_extra_artists=bbox_extra_artists,
        dpi=600,
    )
    plt.close()
