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

    update_database = args.update_database

    classify_structures = args.classify_structures
    cluster_structures = args.cluster_structures

    plot_database = args.plot_database
    pymol_database = args.pymol_database
    table_database = args.table_database

    pdbaa_path = args.pdbaa_path
    nomenclature_path = args.nomenclature_path

    num_cpu = args.num_cpu

    if data_path is not None:
        prep_rascore(data_path=data_path)

    if update_database:
        update_rascore(
            data_path=data_path, pdbaa_fasta_path=pdbaa_path, num_cpu=num_cpu
        )

    if cluster_structures is not None:
        cluster_rascore(
            out_path=out_path,
            data_path=data_path,
            name_table_path=nomenclature_path,
        )

    if classify_structures is not None:
        classify_rascore(
            coord_paths=classify_structures,
            output_path=out_path,
            data_path=data_path,
            num_cpu=num_cpu,
        )

    if plot_database:
        plot_rascore(out_path=out_path, data_path=data_path)
        print("Created rascore database plots!")

    if pymol_database:
        pymol_rascore(out_path=out_path, data_path=data_path)
        print("Created rascore database  PyMOLs!")

    if table_database:
        table_rascore(out_path=out_path, data_path=data_path)
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
        help=f"path to racore database directory (e.g., {os.getcwd()}/{rascore_str}_{data_str})",
    )
    parser.add_argument(
        "-out",
        "--out_path",
        type=str,
        required=False,
        help=f"output path; default is {os.getcwd()}",
    )
    parser.add_argument(
        "-update",
        "--update_database",
        action="store_true",
        required=False,
        help=f"update rascore database from the Protein Data Bank; files saved to {os.getcwd()}/{rascore_str}_{data_str} unless specified otherwise with -data flag",
    )
    parser.add_argument(
        "-classify",
        "--classify_structures",
        nargs="+",
        required=False,
        help=f"space separated list of coordinate file paths (PDB or mmCIF formats) to conformationally classify or path to tab-separated text file with columns (1) core_path (coordinate file path to classify), (2) modelid (0 or another model number to classsify), (3) chainid (chain to classify), and (4) nuc_class (optional to overwride automated nucleotide classification); files saved to {rascore_str}_{classify_str} folder within specified output path",
    )
    parser.add_argument(
        "-cluster",
        "--cluster_structures",
        action="store_true",
        required=False,
        help=f"path to tab-separated text file with columns (1) loop (SW1 or SW2), (2) nucleotide (0P, 2P, or 3P), (3) rama (ABLE sequence), and (4) cluster (conformation name) can be specified with -nomenclature flag; files saved to {rascore_str}_{cluster_str} folder within specified output path",
    )
    parser.add_argument(
        "-plot",
        "--plot_database",
        action="store_true",
        required=False,
        help=f"create plots for rascore database; files saved to {rascore_str}_{plot_str} folder within specified output path",
    )
    parser.add_argument(
        "-pymol",
        "--pymol_database",
        action="store_true",
        required=False,
        help=f"create PyMOL sessions for rascore database; ; files saved to {rascore_str}_{pymol_str} folder within specified output path",
    )
    parser.add_argument(
        "-table",
        "--table_database",
        action="store_true",
        required=False,
        help=f"create tables for rascore database; files saved to {rascore_str}_{table_str} folder within specified output path",
    )
    parser.add_argument(
        "-pdbaa",
        "--pdbaa_path",
        type=str,
        required=False,
        help="path to alternative pdbaa file (for -update flag)",
    )
    parser.add_argument(
        "-nomenclature",
        "--nomenclature_path",
        type=str,
        required=False,
        help="path to alternative nomenclature file (for -cluster flag)",
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