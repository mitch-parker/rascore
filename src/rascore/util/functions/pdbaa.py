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


def get_pdbaa_pdb_code(desc):

    return desc[0:4].lower()


def get_pdbaa_chainid(desc):

    return desc[4:5]


def get_desc_item(desc, index):

    return str(desc).split(None, 8)[index]


def get_pdbaa_method(desc):

    return get_desc_item(desc, 3)


def get_pdbaa_resolution(desc):

    return get_desc_item(desc, 4)


def get_pdbaa_r_factor(desc):

    return get_desc_item(desc, 5)


def get_pdbaa_prot(desc):

    return str(get_desc_item(desc, 8)).split(" <")[0]


def get_pdbaa_swiss_id(desc):

    return str(desc.split("<")[1].split(">")[0].split("(")[0])