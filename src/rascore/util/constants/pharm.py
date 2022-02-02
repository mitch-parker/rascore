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

sp2_name = "SP2"
sp12_name = "SP12"
as_name = "Allosteric Site"
ns_name = "Nucleotide Site"
bs_name = "Base Site"
ps_name = "P110 Site"
other_pharm_name = "Other"
mult_pharm_name = "Multiple"
none_pharm_name = "None"

acr_name = "Acrylamide"
sul_name = "Sulfonamide"

ind_name = "Indole"
ben_name = "Benzodioxane"
bip_name = "Biphenyl"

unclass_name = "Unclassified"

sp2_color = "#d95f02"
sp12_color = "#1b9e77"
other_pharm_color = "#66a61e"
mult_pharm_color = "#e6ab02"
none_pharm_color = "#a6761d"

pharm_class_lst = [
    sp2_name,
    sp12_name,
    other_pharm_name,
    mult_pharm_name,
    none_pharm_name,
]
pharm_color_dict = {
    sp2_name: sp2_color,
    sp12_name: sp12_color,
    other_pharm_name: other_pharm_color,
    mult_pharm_color: mult_pharm_color,
    none_pharm_name: none_pharm_color,
}

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