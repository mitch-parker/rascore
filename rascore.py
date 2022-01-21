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
import argparse
import pyfiglet

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

    print(pyfiglet.figlet_format("rascore"))
    print("Author: Mitchell Parker <mitch.parker@gmail.com>\n")

    if data_path is not None:
        prep_rascore(data_path=data_path)

    if update_database:
        update_rascore(
            data_path=data_path,
            pdbaa_fasta_path=pdbaa_path,
            num_cpu=num_cpu,
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

    classify_str = "three options to input coordinate file paths (PDB or mmCIF formats) for conformational classification"
    classify_str = ";"
    classify_str += f"files saved to {rascore_str}_{classify_str} folder within specified output path"
    classify_str = ":"
    classify_str += "\n"
    classify_str += "a) space separated list "
    classify_str += "\n"
    classify_str += "b) line seperated text file"
    classify_str += "\n"
    classify_str += "c) tab-separated text file with columns (1) core_path (coordinate file path), (2) modelid (optional, 0 or another model number), (3) chainid (chain identifier), and (4) nuc_class (optional, overwrites suggested nucleotide state)"

    parser.add_argument(
        "-classify",
        "--classify_structures",
        nargs="+",
        required=False,
        help=classify_str,
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
        "-edia",
        "--update_edia",
        action="store_true",
        required=False,
        help="update EDIA scores",
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