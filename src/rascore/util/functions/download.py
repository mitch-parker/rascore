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

import urllib.request
import signal

from .path import append_file_path, delete_path, path_exists, unzip_file


def alarm_handler(signum, frame):

    raise Exception


def download_file(url, path, check=True, alarm=60, tries=3):

    append_file_path(path)

    present = False
    if check:
        present = path_exists(path)

    if not present:
        for _ in range(1, tries):
            try:
                signal.signal(signal.SIGALRM, alarm_handler)
                signal.alarm(alarm)
                urllib.request.urlretrieve(url, path)
                signal.alarm(0)
                break
            except Exception:
                delete_path(path)
                signal.alarm(0)
                alarm *= 1.5
            signal.alarm(0)


def download_unzip(url, path, check=True):

    append_file_path(path)

    present = False
    if check:
        present = path_exists(path)

    if not present:
        unzip_file(urllib.request.urlopen(url), out_path=path)
