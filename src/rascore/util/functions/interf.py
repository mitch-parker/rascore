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