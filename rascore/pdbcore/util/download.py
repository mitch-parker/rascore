# -*- coding: utf-8 -*-

"""
Copyright (C) 2021 Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>

This file is part of the rascore project.

The rascore project can not be copied, edited, and/or distributed without the express
permission of Mitchell Isaac Parker <mitch.isaac.parker@gmail.com>.
"""

import urllib.request
import signal

from .path import delete_path, path_exists, append_file_path, unzip_file


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


def download_unzip(url, path):

    append_file_path(path)

    unzip_file(urllib.request.urlopen(url), out_path=path)
