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
