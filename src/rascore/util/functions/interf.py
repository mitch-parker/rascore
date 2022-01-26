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