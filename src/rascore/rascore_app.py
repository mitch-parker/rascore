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

import streamlit as st
from PIL import Image


from scripts import *
from scripts.pages import classify_page, home_page, pdb_page, query_page


class MultiPage:
    def __init__(self) -> None:
        self.pages = []

    def add_page(self, title, func) -> None:
        self.pages.append({"title": title, "function": func})

    def run(self):
        page = st.sidebar.selectbox(
            "Menu", self.pages, format_func=lambda page: page["title"]
        )
        page["function"]()


app = MultiPage()

img = Image.open(
    get_file_path(
        "rascore_logo.png",
        dir_path=get_dir_path(
            dir_str=f"{scripts_str}/{data_str}", dir_path=get_dir_name(__file__)
        ),
    )
)
st.image(img)

app.add_page("Home", home_page.app)
app.add_page("Search for PDB Entry", pdb_page.app)
app.add_page("Query Database", query_page.app)
app.add_page("Classify RAS Structure(s)", classify_page.app)

app.run()
