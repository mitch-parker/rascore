# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from scripts import *
from constants import *


def classify_rascore(coord_path_lst, data_path=None):

    if data_path is None:
        data_path = get_dir_path(dir_str="data", dir_path=get_dir_name(__file__))