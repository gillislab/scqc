#!/usr/bin/env python


import argparse
import io
import itertools
import json
import logging
import os
import gzip
import ast
import glob
from pathlib import Path

# from queue import Queue

import requests
import subprocess
import sys
import time

from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty

import xml.etree.ElementTree as et
import pandas as pd

from scqc.utils import *


# Translate between Python and SRAToolkit log levels for wrapped commands.
#  fatal|sys|int|err|warn|info|debug
LOGLEVELS = {
    10: 'debug',
    20: 'info',
    30: 'warn',
    40: 'err',
    50: 'fatal',
}


class AssignCellType(object):
    def __init__(self, config, soloout_dir,  outlist):
        self.log = logging.getLogger('sra')
        self.config = config
        self.outlist = outlist
        self.max_rank == config.get('metamarker', 'max_rank')

        self.soloout_dir = soloout_dir
        self.outdir = os.path.expanduser(
            self.config.get('metamarker', 'outdir'))
        self.marker_dir = os.path.expanduser(
            self.config.get('metamarker', 'marker_dir'))
        self.sracache = os.path.expanduser(
            self.config.get('metamarker', 'cachedir'))

    def load_marker_sets(self):
        # run R script
        cmd = ["Rscript",
               "bin/getMarkers.R",
               "--marker_direc", f'{self.marker_dir}',
               ]
        subprocess.run(cmd)
        pass

    def assign_cell_types():
        # run R script
        pass

    def execute():
        pass
