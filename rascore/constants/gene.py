# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

hras_name = "HRAS"
kras_name = "KRAS"
nras_name = "NRAS"

swiss_id_lst = ["RASK_HUMAN", "RASN_HUMAN", "RASH_HUMAN"]
uniprot_id_lst = ["P01116", "P01116-2", "P01112", "P01111"]

gene_class_lst = [hras_name, kras_name, nras_name]

gene_class_dict = {
    "GTPase HRas": hras_name,
    "GTPase KRas": kras_name,
    "GTPase NRas": nras_name,
}