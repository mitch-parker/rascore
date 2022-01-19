# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""
from ..functions import *
from .nuc import *
from .pharm import *
from .prot import *

sw1_name = "SW1"
sw2_name = "SW2"

sw1_resids = "25-40"
sw2_resids = "56-76"

loop_resid_dict = {sw1_name: sw1_resids, sw2_name: sw2_resids}

sw1_color = pink_hex
sw2_color = blue_hex

y32_name = "Y32"
y71_name = "Y71"

sw1_nf_name = f"{sw1_name}.{nf_name}"
sw1_gdp_name = f"{sw1_name}.{gdp_name}"
sw1_gtp_name = f"{sw1_name}.{gtp_name}"

sw2_nf_name = f"{sw2_name}.{nf_name}"
sw2_gdp_name = f"{sw2_name}.{gdp_name}"
sw2_gtp_name = f"{sw2_name}.{gtp_name}"

wat_hb_name = "WaterHB"
dir_hb_name = "DirectHB"
no_hb_name = "NoHB"

r_name = "R"
t_name = "T"

noise_name = "Noise"

sw1_gtp_wat_name = f"{sw1_gtp_name}-{wat_hb_name}"
sw1_gtp_dir_name = f"{sw1_gtp_name}-{dir_hb_name}"
sw1_gtp_no_name = f"{sw1_gtp_name}-{no_hb_name}"

sw1_gtp_wat_color = cyan_hex
sw1_gtp_dir_color = purple_hex
sw1_gtp_no_color = pink_hex

noise_color = gray_hex

sw1_gtp_dict = {
    wmhb_name: sw1_gtp_wat_name,
    hb_name: sw1_gtp_dir_name,
    no_hb_name: sw1_gtp_no_name,
}

sw1_gtp_color_dict = {
    sw1_gtp_wat_name: sw1_gtp_wat_color,
    sw1_gtp_dir_name: sw1_gtp_dir_color,
    sw1_gtp_no_name: sw1_gtp_no_color,
}


sw2_nf_gef_name = f"{sw2_nf_name}-{gef_name}"

sw2_gdp_sp12_name = f"{sw2_gdp_name}-{sp12_name}"
sw2_gdp_sp2_a_name = f"{sw2_gdp_name}-{sp2_name}-A"
sw2_gdp_sp2_b_name = f"{sw2_gdp_name}-{sp2_name}-B"
sw2_gdp_binder_name = f"{sw2_gdp_name}-{binder_name}"

sw2_gtp_r_name = f"{sw2_gtp_name}-{r_name}"
sw2_gtp_sp12_a_name = f"{sw2_gtp_name}-{sp12_name}-A"
sw2_gtp_sp12_b_name = f"{sw2_gtp_name}-{sp12_name}-B"
sw2_gtp_t_name = f"{sw2_gtp_name}-{t_name}"

conf_name_dict = {
    sw1_name: {
        nf_name: {"BBAABBBBBAABBLAA": sw1_nf_name},
        gdp_name: {
            "ALBBBABBBAABBBBB": sw1_gdp_name,
            "ALBBBAABBAABBBBB": sw1_gdp_name,
            "ALBBBBLBBAABBBBB": sw1_gdp_name,
        },
        gtp_name: {
            "ALBBBBBBBBBABBBB": sw1_gtp_name,
            "ALBBBABBBBBABBBB": sw1_gtp_name,
            "ALBBBABBBBBALBBB": sw1_gtp_name,
        },
    },
    sw2_name: {
        nf_name: {"BBABLAAABAAAAAAAAAABAm": sw2_nf_gef_name},
        gdp_name: {
            "BBBAEBBABAAAAAAAAAABAp": sw2_gdp_sp12_name,
            "BBBALBBABAAAAAAAAAABAp": sw2_gdp_sp12_name,
            "BBBABBBABAAAAAAAAAABAp": sw2_gdp_sp12_name,
            "BBBAABBABAAAAAAAAAABAp": sw2_gdp_sp12_name,
            "BBBBEBBBBBAAAAAAAAABAm": sw2_gdp_sp2_a_name,
            "BBBBEBBBBBAAAAAAAAABAt": sw2_gdp_sp2_a_name,
            "BBBALBABBBAAAAAAAAABAm": sw2_gdp_sp2_b_name,
            "BBBBLAABBBAAAAAAAAABAm": sw2_gdp_binder_name,
        },
        gtp_name: {
            "BBBBABAAAAAAAAAAAAABAp": sw2_gtp_r_name,
            "BBBBABAABAAAAAAAAAABAp": sw2_gtp_r_name,
            "BBBBABBABAAAAAAAAAABAp": sw2_gtp_sp12_a_name,
            "BBBBAABBBBAAAAAAAAABAp": sw2_gtp_sp12_b_name,
            "BBBBABBBBBAAAAAAAAABAp": sw2_gtp_sp12_b_name,
            "BBBBABBABBAAAAAAAAABAp": sw2_gtp_sp12_b_name,
            "BBBBABBBBBAAAAAABAABAp": sw2_gtp_sp12_b_name,
            "BBBBABAAAAAAAAAAAAABAm": sw2_gtp_t_name,
            "BBBBABAAAAABAAAAAAABAm": sw2_gtp_t_name,
            "BBBBABAAAABAAAAAAAABAm": sw2_gtp_t_name,
            "BBBBABAAAABLAAAAAAABAm": sw2_gtp_t_name,
        },
    },
}


conf_color_dict = {
    sw1_name: {
        sw1_nf_name: nf_color,
        sw1_gdp_name: gdp_color,
        sw1_gtp_wat_name: sw1_gtp_wat_color,
        sw1_gtp_dir_name: sw1_gtp_dir_color,
        sw1_gtp_no_name: sw1_gtp_no_color,
        noise_name: noise_color,
    },
    sw2_name: {
        sw2_nf_gef_name: light_blue_hex,
        sw2_gdp_sp12_name: light_orange_hex,
        sw2_gdp_sp2_a_name: light_olive_hex,
        sw2_gdp_sp2_b_name: light_red_hex,
        sw2_gdp_binder_name: light_purple_hex,
        sw2_gtp_r_name: light_green_hex,
        sw2_gtp_sp12_a_name: light_cyan_hex,
        sw2_gtp_r_name: light_pink_hex,
        sw2_gtp_sp12_b_name: light_brown_hex,
        noise_name: light_gray_hex,
    },
}


conf_nuc_color_dict = {
    sw1_name: {
        sw1_nf_name: nf_color,
        sw1_gdp_name: gdp_color,
        sw1_gtp_name: gtp_color,
        noise_name: noise_color,
    },
    sw2_name: {
        sw2_nf_gef_name: nf_color,
        sw2_gdp_sp12_name: gdp_color,
        sw2_gdp_sp2_a_name: gdp_color,
        sw2_gdp_sp2_b_name: gdp_color,
        sw2_gdp_binder_name: gdp_color,
        sw2_gtp_r_name: gtp_color,
        sw2_gtp_sp12_a_name: gtp_color,
        sw2_gtp_t_name: gtp_color,
        sw2_gtp_sp12_b_name: gtp_color,
        noise_name: noise_color,
    },
}