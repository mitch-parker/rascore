# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import math

from .lst import lst_unique


def calc_q_score(i_cont_lst, j_cont_lst, i_dist_lst, j_dist_lst):

    cont_lst = lst_unique(i_cont_lst, j_cont_lst)

    q_sum = 0
    w_sum = 0

    for cont in cont_lst:

        i_cb_dist = 12
        if cont in i_cont_lst:
            i_cb_dist = i_dist_lst[i_cont_lst.index(cont)]

        j_cb_dist = 12
        if cont in j_cont_lst:
            j_cb_dist = j_dist_lst[j_cont_lst.index(cont)]

        min_cb_diff = min(i_cb_dist, j_cb_dist)
        abs_cb_diff = abs(i_cb_dist - j_cb_dist)

        if min_cb_diff > 5:
            weight = math.exp(-((min_cb_diff - 5) ** 2) / 9.158)
        else:
            weight = 1

        q_sum += weight * math.exp(-0.5 * abs_cb_diff)
        w_sum += weight

    if w_sum == 0:
        q_score = 1
    else:
        q_score = 1 - (q_sum / w_sum)

    return q_score