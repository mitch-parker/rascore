# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from .scripts import *
from .constants import *


def pymol_rascore(data_path=None, out_path=None):

    if data_path is None:
        data_path = f"{os.getcwd()}/{rascore_str}_{data_str}"
    if out_path is None:
        out_path = f"{os.getcwd()}/{rascore_str}_{pymol_str}"

    print("Rascore PyMOLs complete!")