# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from scripts.functions import *

nf_name = "0P"
gdp_name = "2P"
gtp_name = "3P"

nf_color = blue_hex
gdp_color = orange_hex
gtp_color = green_hex

nuc_color_dict = {nf_name: nf_color, gdp_name: gdp_color, gtp_name: gtp_color}

nuc_class_dict = {
    "GNP": gtp_name,
    "GDP": gdp_name,
    "GTP": gtp_name,
    "GSP": gtp_name,
    "GCP": gtp_name,
    "DBG": gtp_name,
    "AGN": gtp_name,
    "9GM": gtp_name,
    "CAG": gtp_name,
    "None": nf_name,
}
