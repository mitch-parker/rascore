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
    #cluster = args.cluster
    gui = args.gui

    out = args.out
    cpu = args.cpu

    if classify is not None:

        from rascore.util.pipelines.classify_rascore import classify_rascore

        classify_rascore(
            file_paths=classify,
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

    # elif cluster is not None:
    #     from util.pipelines.prep_rascore import prep_rascore
    #     from util.pipelines.cluster_rascore import cluster_rascore

    #     prep_rascore(build_path=cluster)
    #     cluster_rascore(cluster, out_path=out, num_cpu=cpu)

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

    print(pyfiglet.figlet_format("Rascore"))
    print("A tool for analyzing RAS protein structures\n")
    print("Author: Mitchell Isaac Parker <mip34@drexel.edu>")
    print("License: Apache License 2.0\n")

    if not args:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Rascore: A tool for analyzing RAS protein structures"
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-classify",
        "--classify",
        nargs="+",
        required=False,
        help="path to coordinate files of RAS structures to conformationally classify provided as a space-separated list, line-sepated file, or tab-seperated file (output files saved to rascore_classify in current working directory unless an output directory path is specified)",
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
    # group.add_argument(
    #     "-cluster",
    #     "--cluster",
    #     type=str,
    #     required=False,
    #     help="path to rascore database directory (output files saved to rascore_cluster in current working directory unless an output directory path is specified)",
    # )
    group.add_argument(
        "-gui",
        "--gui",
        type=str,
        const=True,
        nargs="?",
        required=False,
        help="optional path to rascore database directory (can run without) for running the rascore GUI application",
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
