# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Mitchell Isaac Parker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from ..functions.color import (
    gray_hex,
    pink_hex,
    purple_hex,
    cyan_hex,
    olive_hex,
    light_blue_hex,
    light_orange_hex,
    light_olive_hex,
    light_red_hex,
    light_purple_hex,
    light_green_hex,
    light_cyan_hex,
    light_pink_hex,
    light_brown_hex,
)
from .nuc import nf_name, gdp_name, gtp_name, nf_color, gdp_color, gtp_color, nuc_class_lst
from .pharm import sp2_name, sp12_name
from .prot import gef_name, binder_name

y32_name = "Y32"
q61_name = "Q61"
y71_name = "Y71"

sw1_name = "SW1"
sw2_name = "SW2"

sw1_resids = "25-40"
sw2_resids = "56-76"

sw1_color = "#e7298a"
sw2_color = "#7570b3"

in_name = 'in'
out_name = 'out'

loop_resid_dict = {sw1_name: sw1_resids, sw2_name: sw2_resids}
loop_color_dict = {sw1_name: sw1_color, sw2_name: sw2_color}

on_name = 'On'
off_name = 'Off'

r_name = "R"
t_name = "T"

outlier_name = 'Outlier'
disorder_name = 'Disordered'

noise_name = "Noise"

outlier_color = gray_hex
disorder_color = gray_hex

y32_in_name = f"{y32_name}{in_name}"
y32_out_name = f"{y32_name}{out_name}"
y71_in_name = f"{y71_name}{in_name}"
y71_out_name = f"{y71_name}{out_name}"


y32_in_gtp_name = f"{y32_in_name}.{gtp_name}"
y32_out_gtp_name = f"{y32_out_name}.{gtp_name}"

y32_in_gdp_name = f"{y32_in_name}.{gdp_name}"
y32_out_gdp_name = f"{y32_out_name}.{gdp_name}"

y32_in_nf_name = f"{y32_in_name}.{nf_name}"
y32_out_nf_name = f"{y32_out_name}.{nf_name}"

y71_in_gtp_name = f"{y71_in_name}.{gtp_name}"
y71_out_gtp_name = f"{y71_out_name}.{gtp_name}"

y71_in_gdp_name = f"{y71_in_name}.{gdp_name}"
y71_out_gdp_name = f"{y71_out_name}.{gdp_name}"

y71_in_nf_name = f"{y71_in_name}.{nf_name}"
y71_out_nf_name = f"{y71_out_name}.{nf_name}"

sw1_gtp_in_on_name = f"{y32_in_gtp_name}-{on_name.upper()}"
sw1_gdp_out_off_name = f"{y32_out_gdp_name}-{off_name.upper()}"
sw1_nf_out_gef_name = f"{y32_out_nf_name}-{gef_name}"

sw2_gtp_in_r_name = f"{y71_in_gtp_name}-{r_name}"
sw2_gtp_in_sp12a_name = f"{y71_in_gtp_name}-{sp12_name}-A"
sw2_gtp_in_sp12b_name = f"{y71_in_gtp_name}-{sp12_name}-B"
sw2_gtp_out_t_name = f"{y71_out_gtp_name}-{t_name}"
sw2_gdp_in_sp12_name = f"{y71_in_gdp_name}-{sp12_name}"
sw2_gdp_out_sp2a_name = f"{y71_out_gdp_name}-{sp2_name}-A"
sw2_gdp_out_sp2b_name = f"{y71_out_gdp_name}-{sp2_name}-B"
sw2_gdp_out_binder_name = f"{y71_out_gdp_name}-{binder_name.upper()}"
sw2_nf_out_gef_name = f"{y71_out_nf_name}-{gef_name}"

conf_name_dict = {
    y32_in_gtp_name: {"ALBBBABBBBBABBBB": sw1_gtp_in_on_name},
    y32_out_gdp_name: {"ALBBBABBBAABBBBB": sw1_gdp_out_off_name}, 
    y32_out_nf_name: {"BBAABBBBBAABBLAA": sw1_nf_out_gef_name},
    y71_in_gtp_name: {"BBBBABAAAAAAAAAAAAABA": sw2_gtp_in_r_name,
                    "BBBBABBABAAAAAAAAAABA": sw2_gtp_in_sp12a_name,
                    "BBBBABBBBBAAAAAAAAABA": sw2_gtp_in_sp12b_name},
    y71_out_gtp_name:{"BBBBABAAAABLAAAAAAABA": sw2_gtp_out_t_name},
    y71_in_gdp_name: {"BBBAEBBABAAAAAAAAAABA": sw2_gdp_in_sp12_name},
    y71_out_gdp_name: {"BBBBEBBBABAAAAAAAAABA": sw2_gdp_out_sp2a_name,
                    "BBBALBABBBAAAAAAAAABA": sw2_gdp_out_sp2b_name,
                    "BBBBLAABBBAAAAAAAAABA": sw2_gdp_out_binder_name},
    y71_out_nf_name: {"BBABLAAABAAAAAAAAAABA": sw2_nf_out_gef_name}
}

conf_color_dict = {
    sw1_name: {
        sw1_nf_out_gef_name: nf_color,
        sw1_gdp_out_off_name: gdp_color,
        sw1_gtp_in_on_name: gtp_color,
        disorder_name: disorder_color,
        outlier_name: outlier_color
    },
    sw2_name: {
        sw2_nf_out_gef_name: light_blue_hex,
        sw2_gdp_out_sp2a_name: light_orange_hex,
        sw2_gdp_in_sp12_name: light_olive_hex,
        sw2_gdp_out_sp2b_name: light_red_hex,
        sw2_gdp_out_binder_name: light_purple_hex,
        sw2_gtp_in_r_name: light_green_hex,
        sw2_gtp_in_sp12a_name: light_cyan_hex,
        sw2_gtp_out_t_name: light_pink_hex,
        sw2_gtp_in_sp12b_name: light_brown_hex,
        disorder_name: disorder_color,
        outlier_name: outlier_color
    },
}

conf_nuc_color_dict = {
    sw1_name: {
        sw1_nf_out_gef_name: nf_color,
        sw1_gdp_out_off_name: gdp_color,
        sw1_gtp_in_on_name: gtp_color,
        outlier_name: outlier_color,
        disorder_name: disorder_color,
    },
    sw2_name: {
        sw2_nf_out_gef_name: nf_color,
        sw2_gdp_out_sp2a_name: gdp_color,
        sw2_gdp_in_sp12_name: gdp_color,
        sw2_gdp_out_sp2b_name: gdp_color,
        sw2_gdp_out_binder_name: gdp_color,
        sw2_gtp_in_r_name: gtp_color,
        sw2_gtp_in_sp12a_name: gtp_color,
        sw2_gtp_out_t_name: gtp_color,
        sw2_gtp_in_sp12b_name: gtp_color,
        outlier_name: outlier_color,
        disorder_name: disorder_color,
    },
}


resid_color_dict = {y32_name:{y32_in_name: pink_hex, y32_out_name: cyan_hex},
                y71_name:{y71_in_name: purple_hex, y71_out_name: olive_hex}}