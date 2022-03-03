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

import numpy as np


def format_lst(val_lst, return_str=False, return_int=False, return_float=False):

    return_type = type(val_lst)

    return return_type(
        [
            format_val(
                x,
                return_str=return_str,
                return_int=return_int,
                return_float=return_float,
            )
            for x in val_lst
        ]
    )


def format_val(val, return_str=False, return_int=False, return_float=False):

    if type(val) == list or type(val) == tuple:
        val = format_lst(
            val, return_str=return_str, return_int=return_int, return_float=return_float
        )
    else:
        if return_str and not return_int and not return_float:
            val = str(val)
        elif return_int and not return_str and not return_float:
            try:
                val = int(val)
            except Exception:
                pass
        elif return_float and not return_str and not return_int:
            try:
                val = float(val)
            except Exception:
                pass

    return val


def format_nested_lst(val_lst, return_str=False, return_int=False, return_float=False):

    for i, val in enumerate(val_lst):
        val_lst[i] = [
            format_val(
                x,
                return_str=return_str,
                return_int=return_int,
                return_float=return_float,
            )
            for x in val
        ]

    return val_lst


def type_lst(data, return_str=False, return_int=False, return_float=False, sort=False):

    if type(data) != list:
        data = list([data])

    data = format_lst(
        data, return_str=return_str, return_int=return_int, return_float=return_float
    )

    if sort:
        data = sorted(data)

    return data


def lst_unique(starting_lst, return_str=False, return_int=False, return_float=False):

    val_lst = list()
    [val_lst.append(x) for x in starting_lst if x not in val_lst]

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def move_end_lst(val_lst, end_lst):

    end_lst = type_lst(end_lst)

    for end in end_lst:
        if end in val_lst:
            val_lst.append(val_lst.pop(val_lst.index(end)))

    return val_lst


def sort_lst(val_lst, return_str=False, return_int=False, return_float=False):

    val_dict = dict()
    for index, val in enumerate(val_lst):
        try:
            fix_val = int(val)
        except:
            fix_val = str(val)

        val_lst[index] = fix_val
        val_dict[fix_val] = val

    val_lst = sorted(
        val_lst,
        key=lambda v: (
            isinstance(v, str),
            v,
        ),
    )

    val_lst = [val_dict[x] for x in val_lst]

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def add_lsts(lst_1, lst_2, return_str=False, return_int=False, return_float=False):

    val_lst = list(set(lst_1 + lst_2))

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def subtract_lsts(lst_1, lst_2, return_str=False, return_int=False, return_float=False):

    val_lst = list(set(lst_1) - set(lst_2))

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def lst_inter(lst_1, lst_2, return_str=False, return_int=False, return_float=False):

    val_lst = list(set(lst_1).intersection((set(lst_2))))

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def lst_diff(
    right_lst, left_lst, return_str=False, return_int=False, return_float=False
):

    val_lst = list(set(left_lst).difference((set(right_lst))))

    val_lst = format_lst(
        val_lst, return_str=return_str, return_int=return_int, return_float=return_float
    )

    return val_lst


def calc_jaccard(lst_1, lst_2, return_dist=False):

    inter = len(lst_inter(lst_1, lst_2))
    union = (len(lst_1) + len(lst_2)) - inter

    jaccard = float(inter) / union

    if return_dist:
        jaccard = 1 - jaccard

    return jaccard


def calc_simpson(lst_1, lst_2, return_dist=False):

    intersect = len([x for x in lst_1 if x in lst_2])

    min_size = np.min(np.array([len(lst_1), len(lst_2)]))

    simpson = intersect / min_size

    if return_dist:
        simpson = 1 - simpson

    return simpson


def lst_nums(start, end, return_str=False, return_int=False, return_float=False):

    first = int(start)
    last = int(end) + 1

    num_lst = list(range(first, last))

    num_lst = format_lst(
        num_lst,
        return_str=return_str,
        return_int=return_int,
        return_float=return_float,
    )

    return num_lst


def lst_to_str(val_lst, join_txt=",", empty=None):

    val_lst = type_lst(val_lst)

    if len(val_lst) == 0:
        val_str = str(empty)
    elif len(val_lst) == format_val(val_lst[0], return_str=True) == "None":
        val_str = str(empty)
    else:
        if len(val_lst) == 0:
            val_str = str(empty)
        elif len(val_lst) > 0:

            val_lst = format_lst(val_lst, return_str=True)
            val_str = str(join_txt.join(val_lst))

    return val_str


def str_to_lst(
    val_str,
    sep_txt=",",
    return_str=False,
    return_int=False,
    return_float=False,
):

    val_str = format_val(val_str, return_str=True)

    if sep_txt in val_str:
        val_lst = list(val_str.split(sep_txt))
    else:
        val_lst = type_lst(val_str)

    val_lst = format_lst(
        val_lst,
        return_str=return_str,
        return_int=return_int,
        return_float=return_float,
    )

    return val_lst


def res_to_str(res_lst, sep_txt=":"):

    if res_lst is None:
        res_str = "all"
    else:
        res_lst = type_lst(res_lst, return_int=True, sort=True)
        res_str = ""
        count = 0

        prev_res = None

        for index, curr_res in enumerate(res_lst):

            if prev_res is not None:
                if index == (len(res_lst) - 1):
                    res_str += f"-{curr_res}"
                else:
                    if (curr_res - prev_res) == 1:
                        count += 1
                    else:
                        if count != 0:
                            res_str += f"-{prev_res}"
                        res_str += sep_txt
                        count = 0

            if count == 0 and index != (len(res_lst) - 1):
                res_str += str(curr_res)

            prev_res = curr_res

    return res_str


def res_to_lst(res_str, return_str=False):

    if res_str == None or type(res_str) == list:
        res_lst = res_str
    elif type(res_str) == int:
        res_lst = type_lst(res_str)
    elif type(res_str) == str:
        if ":" in res_str:
            val_lst = res_str.split(":")
        else:
            val_lst = type_lst(res_str)

        res_lst = list()

        for val in val_lst:
            if "-" in val:
                res_range = val.split("-")
                res_start = res_range[0]
                res_end = res_range[1]
                res_lst += lst_nums(res_start, res_end)
            else:
                res_lst.append(int(val))

    return_int = True
    if return_str:
        return_int = False

    if res_lst is not None:
        res_lst = format_lst(res_lst, return_str=return_str, return_int=return_int)

    return res_lst


def build_range_lst(range, step, sep="-", type=float, dec=2):

    range_lst = list()
    if sep in str(range):
        sep_range = range.split(sep)
        start = type(sep_range[0])
        end = type(sep_range[1])
        while start <= end:
            range_lst.append(start)
            start += type(step)
        if end not in range_lst:
            range_lst.append(end)
    else:
        range_lst.append(type(range))

    if type == float:
        range_lst = [round(x, dec) for x in range_lst]

    return lst_unique(range_lst)


def get_lst_val_indices(lst, val):

    i = -1
    loc_lst = []
    while True:
        try:
            loc = lst.index(val, i + 1)
        except ValueError:
            break
        else:
            loc_lst.append(loc)
            i = loc
    return loc_lst
