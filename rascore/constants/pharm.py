# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

sp2_name = "SP2"
sp12_name = "SP12"
np_name = "NP"
bp_name = "BP"
p110_name = "P110"
other_pharm_name = "Other"
none_pharm_name = "None"

acr_name = "Acrylamide"
sul_name = "Sulfonamide"

ind_name = "Indole"
ben_name = "Benzodioxane"
bip_name = "Biphenyl"

unclass_name = "Unclassified"

sp2_color = "#d95f02"
sp12_color = "#1b9e77"

pharm_class_lst = [sp2_name, sp12_name, other_pharm_name, none_pharm_name]
pharm_color_dict = {sp2_name: sp2_color, sp12_name: sp12_color}

sp2_cont = [12, 96, 99]
sp12_cont = [5, 39, 54]
np_cont = [29, 30, 32]
bp_cont = [85, 118, 119]
p110_cont = [106, 108, 110]

pharm_site_dict = {
    sp2_name: sp2_cont,
    sp12_name: sp12_cont,
    np_name: np_cont,
    bp_name: bp_cont,
    p110_name: p110_cont,
}

pharm_match_dict = {
    sp2_name: {
        acr_name: ["CC(=O)N1CCNCC1", "CC(=O)N1CCCC1", "CC(=O)N1CCC1"],
        sul_name: ["CCS(=O)(=O)N", "CC(=O)N1CCCCC1"],
    },
    sp12_name: {
        ind_name: ["C1=CC=C2C(=C1)C=CN2"],
        ben_name: ["C1COC2=CC=CC=C2O1"],
        bip_name: ["C1=CC=C(C=C1)C2=CC=CC=C2"],
    },
}