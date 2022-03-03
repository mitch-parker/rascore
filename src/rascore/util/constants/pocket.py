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
from ..scripts.prep_pocket import pocket_unbound_name
from .pharm import sp12_name, sp2_name, other_pharm_name

pocket_site_lst = [sp12_name, sp2_name, other_pharm_name]

pocket_color_dict = {
    sp12_name: "#1b9e77",
    sp2_name: "#d95f02",
    pocket_unbound_name: gray_hex,
}