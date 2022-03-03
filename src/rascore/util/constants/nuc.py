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

from ..functions.color import blue_hex, orange_hex, green_hex

nf_name = "0P"
gdp_name = "2P"
gtp_name = "3P"

rep_dict = {
    "4q21A": gdp_name,
    "1bkdR": nf_name,
    "5p21A": gtp_name,
}

nf_color = blue_hex
gdp_color = orange_hex
gtp_color = green_hex

nuc_class_lst = [nf_name, gdp_name, gtp_name]

nuc_color_dict = {nf_name: nf_color, gdp_name: gdp_color, gtp_name: gtp_color}

nuc_class_dict = {
    "GDP": gdp_name,
    "None": nf_name,
}

gtp_atomids = ["O1G", "O2G", "O3G", "S1G", "S2G", "S3G"]
