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

from ..functions.color import gray_hex

sp2_name = "SP2"
sp12_name = "SP12"
as_name = "Allosteric Site"
ns_name = "Nucleotide Site"
bs_name = "Base Site"
ps_name = "P110 Site"
other_pharm_name = "Other"
mult_pharm_name = "Multiple"
none_pharm_name = "None"

bound_name = "Bound"
unbound_name = "Unbound"

acr_name = "Acrylamide"
sul_name = "Sulfonamide"

ind_name = "Indole"
ben_name = "Benzodioxane"
bip_name = "Biphenyl"

uncl_name = "Unclassified"

sp2_color = "#d95f02"
sp12_color = "#1b9e77"
other_pharm_color = "#66a61e"
mult_pharm_color = "#e6ab02"
none_pharm_color = gray_hex

pharm_color_dict = {
    sp2_name: sp2_color,
    sp12_name: sp12_color,
    mult_pharm_name: mult_pharm_color,
    other_pharm_name: other_pharm_color,
    none_pharm_name: none_pharm_color
}

pharm_class_lst = [
    sp12_name,
    sp2_name,
    other_pharm_name,
    mult_pharm_name,
    none_pharm_name,
]

match_class_lst = [   
    f"{sp12_name}.{ind_name}",
    f"{sp12_name}.{ben_name}",
    f"{sp12_name}.{bip_name}",
    f"{sp12_name}.{uncl_name}",
    f"{sp2_name}.{acr_name}",
    f"{sp2_name}.{sul_name}",
    f"{sp2_name}.{uncl_name}",
    other_pharm_name,
    mult_pharm_name,
    none_pharm_name
]

pocket_class_lst = list()

for pharm_class in pharm_class_lst:
    pocket_class_lst.append(f"{pharm_class}-{bound_name}")
    pocket_class_lst.append(f"{pharm_class}-{unbound_name}")

sp2_cont = [12, 96, 99]
sp12_cont = [5, 39, 54]
ns_cont = [29, 30, 32]
bs_cont = [85, 118, 119]
as_cont = [4, 49, 164]
ps_cont = [106, 108, 110]

pharm_site_dict = {
    sp2_name: sp2_cont,
    sp12_name: sp12_cont,
    ns_name: ns_cont,
    bs_name: bs_cont,
    as_name: as_cont,
    ps_name: ps_cont,
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
