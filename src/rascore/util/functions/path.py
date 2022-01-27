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
import shutil
import gzip
import json
import pandas as pd
import numpy as np

from .table import order_cols, order_rows, get_col_order

rascore_str = "rascore"
build_str = "build"
classify_str = "classify"
cluster_str = "cluster"

util_str = "util"
pipelines_str = "pipelines"
data_str = "data"
pages_str = "pages"
functions_str = "functions"

pdbaa_str = "pdbaa"
core_str = "core"
rcsb_str = "rcsb"
renum_str = "renum"
rcsb_assembly_str = "rcsb_assembly"
renum_assembly_str = "renum_assembly"
sifts_str = "sifts"
edia_str = "edia"
seq_str = "sequence"
interf_str = "interface"
pocket_str = "pocket"
lig_str = "ligand"
eds_str = "eds"
map_str = "2mFo-DFc"
diff_str = "mFo-DFc"
sup_str = "super"


def path_exists(path):

    exists = False
    if path is not None:
        if os.path.isfile(path):
            exists = True
        elif os.path.isdir(path):
            exists = True

    return exists


def append_path(path):

    if not path_exists(path):
        os.makedirs(path)


def get_dir_name(dir_path):

    if "/" in dir_path:
        dir_name = dir_path.rsplit("/", 1)[0]
    else:
        dir_name = os.getcwd()

    return dir_name


def get_file_name(path):

    if "/" in path:
        file_name = path.rsplit("/", 1)[1]
    else:
        file_name = path

    return file_name


def append_file_path(path):

    append_path(get_dir_name(path))


def delete_path(path):

    if path_exists(path):
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def copy_path(source_path, dest_path):

    if path_exists(dest_path):
        delete_path(dest_path)
    shutil.copyfile(source_path, dest_path)


def save_table(path, df, sep="\t", header=True, index=False, fillna="None"):

    append_file_path(path)

    df = df.fillna(fillna)

    df = order_cols(df, get_col_order(df))

    df = order_rows(df)

    df.to_csv(path, sep=sep, header=header, index=index)


def load_table(path, sep="\t"):

    if path_exists(path):
        df = pd.read_csv(path, sep=sep, dtype=str)
        df = order_rows(df)
    else:
        df = None

    return df


def save_matrix(path, matrix, delim=","):

    append_file_path(path)

    np.savetxt(path, matrix, delimiter=delim)


def load_matrix(path, delim=","):

    if path_exists(path):
        matrix = np.loadtxt(path, delimiter=delim)
    else:
        matrix = None

    return matrix


def save_lst(path, val_lst):

    append_file_path(path)

    with open(path, "w") as file:
        for val in val_lst:
            file.write(f"{val}\n")


def load_lst(path):

    if path_exists(path):
        with open(path, "r") as file:
            line_lst = file.read().splitlines()
    else:
        line_lst = None

    return line_lst


def save_json(path, json_dict):

    append_file_path(path)

    with open(path, "w") as file:
        json.dump(json_dict, file)


def load_json(path):

    if path_exists(path):
        with open(path, "r") as file:
            json_dict = json.load(file)
    else:
        json_dict = None

    return json_dict


def unzip_file(in_path, out_path=None):

    if out_path is None:
        out_path = in_path.replace(".gz", "")

    with gzip.open(in_path, "rb") as file_in:
        with open(out_path, "wb") as file_out:
            shutil.copyfileobj(file_in, file_out)


def search_dir(dir_path, file_str):

    return [x for x in os.listdir(dir_path) if file_str in x]


def get_dir_path(dir_str=None, dir_path=None):

    if dir_path is None:
        dir_path = os.getcwd()

    if dir_str is not None:
        dir_path += f"/{dir_str}"

    return dir_path


def get_file_path(file_name, dir_str=None, dir_path=None, pre_str=True):

    file_path = get_dir_path(dir_str=dir_str, dir_path=dir_path)
    file_path += "/"
    if pre_str and dir_str != None:
        file_path += dir_str
        file_path += "_"
    file_path += file_name

    return file_path


def modify_coord_path(path, return_pdb=False, add_h=False):

    if return_pdb:
        path = path.replace(".cif", ".pdb")
        if add_h:
            path = path.replace(".pdb", ".h.pdb")

    return path


def get_core_path(
    pdb_code, chainid, modelid=None, dir_path=None, return_pdb=False, add_h=False
):
    pdb_id = f"{pdb_code}{chainid}"

    if modelid is not None:
        pdb_id += str(modelid)

    return modify_coord_path(
        get_file_path(
            f"{pdb_id}_{core_str}.cif",
            dir_str=core_str,
            dir_path=dir_path,
            pre_str=False,
        ),
        return_pdb=return_pdb,
        add_h=add_h,
    )


def get_rcsb_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}.cif.gz", dir_str=rcsb_str, dir_path=dir_path, pre_str=False
    )


def get_sifts_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}.xml.gz", dir_str=sifts_str, dir_path=dir_path, pre_str=False
    )


def get_renum_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}_{renum_str}.cif",
        dir_str=renum_str,
        dir_path=dir_path,
        pre_str=False,
    )


def get_seq_path(uniprot_acc, dir_path=None):

    return get_file_path(
        f"{uniprot_acc}.fasta", dir_str=seq_str, dir_path=dir_path, pre_str=False
    )


def get_edia_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}_{edia_str}.csv", dir_str=edia_str, dir_path=dir_path, pre_str=False
    )


def get_lig_path(lig, dir_path=None):

    return get_file_path(
        f"{lig}_{lig_str}.sdf", dir_str=lig_str, dir_path=dir_path, pre_str=False
    )


def get_eds_map_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}_{map_str}.ccp", dir_str=map_str, dir_path=dir_path, pre_str=False
    )


def get_eds_diff_path(pdb_code, dir_path=None):

    return get_file_path(
        f"{pdb_code}_{diff_str}.ccp4",
        dir_str=diff_str,
        dir_path=dir_path,
        pre_str=False,
    )


def get_interf_path(pdb_code, chainid, interf, dir_path=None, return_pdb=False):

    return modify_coord_path(
        get_file_path(
            f"{pdb_code}_{chainid}_{interf}_{interf_str}.cif",
            dir_str=interf_str,
            dir_path=dir_path,
            pre_str=False,
        ),
        return_pdb=return_pdb,
    )


def get_pocket_path(pdb_id, pocket, dir_path=None):

    return get_file_path(
        f"{pdb_id}_{pocket}_{pocket_str}.pdb",
        dir_str=pocket_str,
        dir_path=dir_path,
        pre_str=False,
    )


def get_neighbor_path(file_path, dir_str, neighbor_str):

    dir_path = get_dir_name(file_path)
    dir_path = dir_path.split(dir_str)[0]
    dir_path += neighbor_str

    return dir_path
