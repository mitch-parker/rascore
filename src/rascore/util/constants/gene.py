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

hras_name = "HRAS"
kras_name = "KRAS"
nras_name = "NRAS"

gene_class_lst = [kras_name, hras_name, nras_name]

swiss_id_lst = ["RASK_HUMAN", "RASN_HUMAN", "RASH_HUMAN"]
uniprot_acc_lst = ["P01116", "P01116-2", "P01112", "P01111"]

gene_class_dict = {
    "GTPase HRas": hras_name,
    "GTPase KRas": kras_name,
    "GTPase NRas": nras_name,
}