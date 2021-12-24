# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the pdbcore project.

The pdbcore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""


def get_pdbaa_pdb_code(desc):

    return desc[0:4].lower()


def get_pdbaa_chainid(desc):

    return desc[4:5]


def get_desc_item(desc, index):

    return str(desc).split(None, 8)[index]


def get_pdbaa_method(desc):

    return get_desc_item(desc, 3)


def get_pdbaa_resolution(desc):

    return get_desc_item(desc, 4)


def get_pdbaa_r_factor(desc):

    return get_desc_item(desc, 5)


def get_pdbaa_prot(desc):

    return str(get_desc_item(desc, 8)).split(" <")[0]


def get_pdbaa_swiss_id(desc):

    return str(desc.split("<")[1].split(">")[0].split("(")[0])