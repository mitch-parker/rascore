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

import os
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import argparse
import sys
import pyfiglet

from rascore.util.functions.path import (
    get_file_path,
    get_dir_name,
    copy_path,
    util_str,
    data_str,
)
from rascore.util.functions.file import entry_table_file


def main(args):

    classify = args.classify
    build = args.build
    gui = args.gui

    out = args.out
    cpu = args.cpu

    if classify is not None:

        from rascore.util.pipelines.classify_rascore import classify_rascore

        classify_rascore(
            coord_paths=classify,
            out_path=out,
            num_cpu=cpu,
        )
    elif build is not None:

        from rascore.util.pipelines.prep_rascore import prep_rascore
        from rascore.util.pipelines.build_rascore import build_rascore

        pdbaa_fasta_path = None
        if build is not True:
            pdbaa_fasta_path = build
        prep_rascore(build_path=out)

        build_rascore(
            out_path=out,
            pdbaa_fasta_path=pdbaa_fasta_path,
            num_cpu=cpu,
        )
    elif gui is not None:
        if gui is not True:

            from rascore.util.pipelines.prep_rascore import prep_rascore

            prep_rascore(build_path=gui)
            entry_table_path = get_file_path(entry_table_file, dir_path=gui)
            copy_path(
                entry_table_path,
                get_file_path(
                    entry_table_file,
                    dir_path=f"{get_dir_name(__file__)}/{util_str}/{data_str}",
                ),
            )

        rascore_app_path = get_file_path(
            "rascore_gui.py",
            dir_path=get_dir_name(__file__),
        )

        os.system(f"streamlit run {rascore_app_path}")


def cli(args=None):

    print(pyfiglet.figlet_format("rascore"))
    print("A tool for analyzing the conformations of RAS structures\n")
    print("Author: Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>")
    print("License: MIT License\n")

    if not args:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="rascore: A tool for analyzing the conformations of RAS structures"
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-classify",
        "--classify",
        nargs="+",
        required=False,
        help="path to coordinate file(s) of RAS structures to conformationally classify provided as a space-separated list, line-sepated file, or tab-seperated file (output files saved to rascore_classify in current working directory unless an output directory path is specified)",
    )
    group.add_argument(
        "-build",
        "--build",
        type=str,
        const=True,
        nargs="?",
        required=False,
        help="build or update rascore database from the Protein Data Bank (output files saved to rascore_build in current working directory unless an output directory path is specified)",
    )
    group.add_argument(
        "-gui",
        "--gui",
        type=str,
        const=True,
        nargs="?",
        required=False,
        help="path to rascore database directory (can run limited version if not specified) for running the rascore GUI application (output files saved to rascore_app in current working directory unless an output directory path is specified)",
    )
    parser.add_argument(
        "-out",
        "--out",
        type=str,
        required=False,
        help=f"output directory path (default = {os.getcwd()}",
    )
    parser.add_argument(
        "-cpu",
        "--cpu",
        type=int,
        default=1,
        required=False,
        help="number of CPUs to use (default = 1)",
    )
    args = parser.parse_args(args)
    main(args)


if __name__ == "__main__":
    cli()
