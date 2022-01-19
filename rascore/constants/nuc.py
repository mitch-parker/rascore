# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from ..functions import *

nf_name = "0P"
gdp_name = "2P"
gtp_name = "3P"

rep_dict = {
    "4q21A": gdp_name,
    "1bkdR": nf_name,
    "5p21A": gtp_name,
}

nf_color = blue_hex
gdp_color = orange_hex
gtp_color = green_hex

nuc_class_lst = [nf_name, gdp_name, gtp_name]

nuc_color_dict = {nf_name: nf_color, gdp_name: gdp_color, gtp_name: gtp_color}

nuc_class_dict = {
    "GDP": gdp_name,
    "None": nf_name,
}
