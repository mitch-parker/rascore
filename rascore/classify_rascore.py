# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from .scripts import *
from .constants import *


def classify_rascore(coord_paths, data_path=None, out_path=None, num_cpu=1):

    if data_path is None:
        data_path = get_dir_path(dir_str=data_str, dir_path=get_dir_name(__file__))

    print("Rascore classification complete!")