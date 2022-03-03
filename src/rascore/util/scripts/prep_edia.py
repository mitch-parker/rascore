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

import time
import pandas as pd
import requests
import json
from tqdm import tqdm
import concurrent.futures

from ..functions.table import merge_dicts
from ..functions.lst import type_lst
from ..functions.download import download_file
from ..functions.path import (
    path_exists,
    save_table,
    get_edia_path,
    append_path,
    get_dir_path,
    save_json,
    edia_str,
)
from ..functions.col import edia_col, b_factor_col


def download_edia_scores(path, pdb_code):

    if not path_exists(path):

        headers = {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }

        data = '{"edia": {"pdbCode": "%s"}}' % pdb_code

        response = requests.post(
            "https://proteins.plus/api/edia_rest", headers=headers, data=data
        )

        wait = 1
        while True:
            try:
                data = json.loads(response.text)
                status = data["status_code"]
                break
            except Exception:
                time.sleep(wait)
                wait *= 1.5

        if status != 400:

            url = data["location"]

            wait = 1
            while True:
                response = requests.get(url)
                data = json.loads(response.text)
                status = data["status_code"]

                if status == 200:
                    url = data["atom_scores"]
                    download_file(url, path)
                    break
                elif status == 400:
                    message = data["message"]
                    if message == "No electron density file available.":
                        df = pd.DataFrame(
                            columns=[
                                "Structure specifier",
                                "Atom name",
                                "Infile id",
                                "Substructure name",
                                "Substructure id",
                                "Chain",
                                "Element",
                                "EDIA",
                                "EDIA fault analysis",
                                "B factor",
                                "Occupancy",
                            ]
                        )
                        save_table(path, df, sep=",", header=True, index=False)
                        break
                elif status == 202 or status == 429:
                    time.sleep(wait)
                    wait *= 1.5


def get_pdb_edia_dict(pdb_code, edia_dir=None, sifts_dict=None):

    edia_path = get_edia_path(pdb_code, dir_path=edia_dir)

    download_edia_scores(edia_path, pdb_code)

    df = pd.read_csv(edia_path, header=0)

    edia_dict = {pdb_code: dict()}

    if len(df) != 0:

        for index in list(df.index.values):

            add_val = True

            chainid = df.at[index, "Chain"]
            resnum = str(df.at[index, "Substructure id"])
            atomid = df.at[index, "Atom name"]
            edia = df.at[index, "EDIA"]
            b_factor = df.at[index, "B factor"]

            if chainid not in list(edia_dict[pdb_code].keys()):
                edia_dict[pdb_code][chainid] = dict()

            if sifts_dict is not None:
                if resnum in list(sifts_dict[pdb_code][chainid].keys()):
                    resnum = sifts_dict[pdb_code][chainid][str(resnum)]

                else:
                    add_val = False

            if add_val:

                if resnum not in list(edia_dict[pdb_code][chainid].keys()):

                    edia_dict[pdb_code][chainid][resnum] = dict()

                edia_dict[pdb_code][chainid][resnum][atomid] = {
                    edia_col: edia,
                    b_factor_col: b_factor,
                }

    return edia_dict


def prep_edia(
    pdb_codes,
    edia_dir=None,
    edia_json_path=None,
    sifts_dict=None,
    num_cpu=1,
):

    append_path(get_dir_path(dir_str=edia_str, dir_path=edia_dir))

    pdb_code_lst = type_lst(pdb_codes)

    edia_dict = dict()

    if num_cpu == 1:
        for pdb_code in tqdm(
            pdb_code_lst,
            desc="Preparing EDIA scores",
            position=0,
            leave=True,
        ):

            edia_dict = merge_dicts(
                [
                    edia_dict,
                    get_pdb_edia_dict(
                        pdb_code, edia_dir=edia_dir, sifts_dict=sifts_dict
                    ),
                ]
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu) as executor:
            job_lst = [
                executor.submit(
                    get_pdb_edia_dict,
                    pdb_code,
                    edia_dir=edia_dir,
                    sifts_dict=sifts_dict,
                )
                for pdb_code in pdb_code_lst
            ]

            for job in tqdm(
                concurrent.futures.as_completed(job_lst),
                desc="Prepared EDIA scores",
                total=len(job_lst),
                miniters=1,
                position=0,
                leave=True,
            ):

                edia_dict = merge_dicts([edia_dict, job.result()])

    print("Prepared EDIA scores!")

    if edia_json_path is not None:
        save_json(edia_json_path, edia_dict)
    else:
        return edia_dict
