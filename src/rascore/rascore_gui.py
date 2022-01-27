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

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import streamlit as st
from PIL import Image

from rascore.util.functions.path import get_dir_name, get_file_path, util_str, data_str
from rascore.util.pages.home_page import home_page
from rascore.util.pages.pdb_page import pdb_page
from rascore.util.pages.query_page import query_page
from rascore.util.pages.classify_page import classify_page


class MultiApp:
    def __init__(self):
        self.apps = []

    def add_app(self, title, func):
        self.apps.append({"title": title, "function": func})

    def run(self):
        st.sidebar.markdown("## Main Menu")
        app = st.sidebar.radio(
            "Select a Page", self.apps, format_func=lambda app: app["title"]
        )
        st.sidebar.markdown("---")
        app["function"]()


app = MultiApp()

img = Image.open(
    get_file_path(
        "rascore_logo.png", dir_path=f"{get_dir_name(__file__)}/{util_str}/{data_str}"
    )
)
st.image(img)

app.add_app("Home Page", home_page)
app.add_app("Search PDB", pdb_page)
app.add_app("Query Database", query_page)
app.add_app("Classify Structures", classify_page)

app.run()
