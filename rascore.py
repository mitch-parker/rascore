# -*- coding: utf-8 -*-

"""
Copyright (C) 2022 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project cannot be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import os
import argparse

from rascore import *


def main(args):

    data_path = args.data_path
    out_path = args.out_path

    classify_structures = args.classify_structures

    update_database = args.update_database
    plot_database = args.plot_database
    pymol_database = args.pymol_database
    table_database = args.table_database

    cluster_update = args.cluster_update

    num_cpu = args.num_cpu

    if data_path is not None:
        prep_rascore(data_path=data_path)

    if update_database:
        update_rascore(data_path=data_path, num_cpu=num_cpu)

    if cluster_update is not None:
        cluster_rascore(
            data_path,
            out_path=out_path,
            name_table_path=cluster_update,
            num_cpu=num_cpu,
        )

    if classify_structures is not None:
        classify_rascore(
            coord_paths=classify_structures,
            output_path=out_path,
            data_path=data_path,
            num_cpu=num_cpu,
        )

    if plot_database:
        plot_rascore(out_path=out_path, data_path=data_path, num_cpu=num_cpu)
        print("Created rascore database plots!")

    if pymol_database:
        pymol_rascore(out_path=out_path, data_path=data_path, num_cpu=num_cpu)
        print("Created rascore database  PyMOLs!")

    if table_database:
        table_rascore(out_path=out_path, data_path=data_path, num_cpu=num_cpu)
        print("Created rascore database tables!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="rascore: A tool for analyzing the conformations of RAS structures"
    )
    parser.add_argument(
        "-data",
        "--data_path",
        type=str,
        required=False,
        help=f"path to racore database directory (i.e., /path/to/rascore_data)",
    )
    parser.add_argument(
        "-out",
        "--out_path",
        type=str,
        default=os.getcwd(),
        required=False,
        help=f"output path (default = {os.getcwd()})",
    )
    parser.add_argument(
        "-classify",
        "--classify_structures",
        nargs="+",
        required=False,
        help="space separated list of coordinate file paths to conformationally classify (PDB or mmCIF formats) ",
    )
    parser.add_argument(
        "-update",
        "--update_database",
        action="store_true",
        required=False,
        help="update rascore database from the Protein Data Bank",
    )
    parser.add_argument(
        "-plot",
        "--plot_database",
        action="store_true",
        required=False,
        help="create plots for rascore database",
    )
    parser.add_argument(
        "-pymol",
        "--pymol_database",
        action="store_true",
        required=False,
        help="create PyMOL sessions for rascore database",
    )
    parser.add_argument(
        "-table",
        "--table_database",
        action="store_true",
        required=False,
        help="create tables for rascore database",
    )
    parser.add_argument(
        "-cluster",
        "--cluster_update",
        required=False,
        help=f"path to updated source cluster file (tab-separated text file with columns: loop, nucleotide, rama, and cluster; input none to get editable text file called nomenclature.txt)",
    )
    parser.add_argument(
        "-cpu",
        "--num_cpu",
        type=int,
        default=1,
        required=False,
        help="number of CPUs to use (default = 1)",
    )
    args = parser.parse_args()
    main(args)