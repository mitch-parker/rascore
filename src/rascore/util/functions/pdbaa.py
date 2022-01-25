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