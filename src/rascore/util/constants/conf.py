# -*- coding: utf-8 -*-
"""
  Copyright 2022 Mitchell Isaac Parker

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

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
from ..functions.coord import hb_name, no_hb_name, wmhb_name

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
        outlier_name: outlier_color,
        disorder_name: disorder_color,
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
        outlier_name: outlier_color,
        disorder_name: disorder_color,
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


resid_color_dict = {y32_name:{y32_in_name: pink_hex, y32_out_name: cyan_hex, disorder_name: disorder_color},
                y71_name:{y71_in_name: purple_hex, y71_out_name: olive_hex, disorder_name: disorder_color}}

hb_color_dict = {wmhb_name: cyan_hex, hb_name: purple_hex, no_hb_name: pink_hex}

y32_name_lst = list(resid_color_dict[y32_name].keys())
y71_name_lst = list(resid_color_dict[y71_name].keys())

sw1_name_lst = list(conf_color_dict[sw1_name].keys())
sw2_name_lst = list(conf_color_dict[sw2_name].keys())