# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

from scripts.functions import *

sp2_name = "SP2"
sp12_name = "SP12"

acr_name = "Acrylamide"
sul_name = "Sulfonamide"

ind_name = "Indole"
ben_name = "Benzodioxane"
bip_name = "Biphenyl"

unclass_name = "Unclassified"

sp2_acr_name = f"{sp2_name}.{acr_name}"
sp2_sul_name = f"{sp2_name}.{sul_name}"
sp2_unclass_name = f"{sp2_name}.{unclass_name}"

sp12_ind_name = f"{sp12_name}.{ind_name}"
sp12_ben_name = f"{sp12_name}.{ben_name}"
sp12_bip_name = f"{sp12_name}.{bip_name}"
sp12_unclass_name = f"{sp12_name}.{unclass_name}"

other_pharm_name = "Other"
none_pharm_name = "None"

pharm_order_lst = [
    sp2_acr_name,
    sp2_sul_name,
    sp2_unclass_name,
    sp12_ind_name,
    sp12_ben_name,
    sp12_bip_name,
    sp12_unclass_name,
    other_pharm_name,
    none_pharm_name,
]

sp2_color = "salmon"
sp12_color = "slateblue"

pharm_color_dict = {sp2_name: sp2_color, sp12_name: sp12_color}

sp2_cont = [12, 96, 99]
sp12_cont = [5, 39, 54]

pharm_site_dict = {
    sp2_name: sp2_cont,
    sp12_name: sp12_cont,
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