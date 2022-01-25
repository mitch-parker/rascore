import os
import sys
from os import listdir
from os.path import isfile, join
from pathlib import Path
import Bio
from Bio.PDB import *
import pandas as pd  # a great dataframe
import numpy as np  # a large collection of high-level mathematical functions to operate on these arrays
import math
import re  # a regular expression operations
import argparse  # parser for arguments
import time  # time module
from datetime import date
import gzip
import shutil  # compressing module
import xml.etree.ElementTree as ET  # xml parser
import urllib  # urllib
from string import punctuation
import requests
import tqdm  # a progress bar
import socket

from functools import partial
from lxml import html
import requests
from concurrent.futures import as_completed, ProcessPoolExecutor, ThreadPoolExecutor

# exception_AccessionIDs = ["P42212", "Q17104", "Q27903", "Q93125", "P03069", "D3DLN9", "Q96UT3", "P0ABE7", "P00192",
#                           "P76805", "Q8XCE3", "P00720", "Q38170", "Q94N07", "P0AEX9", "P02928", "Q2M6S0"]
socket.setdefaulttimeout(660)
current_directory = os.getcwd()


# # maybe unnecessary
# from multiprocessing import Process
# from multiprocessing import Pool
# from multiprocessing.pool import ThreadPool

# current_directory = os.getcwd()
# default_input_path_to_mmCIF = current_directory + "/mmCIF"
# default_input_path_to_PDB = current_directory + "/PDB"
# default_input_path_to_SIFTS = current_directory + "/SIFTS"
# default_output_path_to_mmCIF = current_directory + "/output_mmCIF"
# default_output_path_to_PDB = current_directory + "/output_PDB"
#
# default_input_path_to_mmCIF_assemblies = current_directory + "/mmCIF_assembly"
# default_input_path_to_PDB_assemblies = current_directory + "/PDB_assembly"
# default_output_path_to_mmCIF_assemblies = current_directory + "/output_mmCIF_assembly"
# default_output_path_to_PDB_assemblies = current_directory + "/output_PDB_assembly"

# current_directory = os.getcwd()
# default_mmCIF_num = 50000
# default_PDB_num = 5000
# gzip_mode = "on"
