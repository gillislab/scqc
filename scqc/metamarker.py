#!/usr/bin/env python


import argparse
import io
import logging
import os
from pathlib import Path

# from queue import Queue

import subprocess
import sys

from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty

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


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()


class Worker(Thread):
    '''
    '''

    def __init__(self, q):
        self.q = q
        super(Worker, self).__init__()

    def run(self):
        # Race condition, just try!
        while True:
            try:
                job = self.q.get_nowait()
                job.execute()
                self.q.task_done()
            except Empty:
                return


class SetUp(object):
    def __init__(self, config):
        self.log = logging.getLogger('sra')
        self.config = config

        self.marker_dir = os.path.expanduser(
            self.config.get('metamarker', 'marker_dir'))
        self.rds_path = os.path.expanduser(
            self.config.get('metamarker', 'rds_path'))
        self.bindir = os.path.expanduser(
            self.config.get('metamarker', 'bindir'))

    def execute(self):
        try:
            os.makedirs(self.marker_dir)
        except FileExistsError:
            pass

        scriptpath = f'{self.bindir}/get_marker.R'
        cmd = ["Rscript", f'{scriptpath}',
               "--mode", "setup",
               "--rds_path", f'{self.rds_path}',
               "--markerdir", f'{self.marker_dir}'
               ]

        cmdstr = " ".join(cmd)
        cp = subprocess.run(cmd)


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
        self.bindir = os.path.expanduser(
            self.config.get('metamarker', 'bindir'))

    def execute(self):

        # run R script
        scriptpath = f'{self.bindir}/get_marker.R'
        outprefix = self.soloout_dir.split("_Solo.out")[0]
        cmd = ['Rscript', f'{scriptpath}',
               '--marker_direc', f'{self.marker_dir}',
               '--solo_out_dir', f'{self.soloout_dir}',
               '--max_rank', f'{self.max_rank}',
               '--outprefix', f'{outprefix}'
               ]
        subprocess.run(cmd)

        pass


# to do: include drivers for 10x and ss alignments
if __name__ == "__main__":

    gitpath = os.path.expanduser("~/git/scqc")
    sys.path.append(gitpath)

    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug',
                        action="store_true",
                        dest='debug',
                        help='debug logging')

    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        dest='verbose',
                        help='verbose logging')

    parser.add_argument('-s', '--setup',
                        action='store_true',
                        dest='setup',
                        help='Set up directories and downloads supplemental data')

    parser.add_argument('-a', '--assign',
                        action='store_true',
                        dest='assign',
                        help='Assign cell types to STAR outout')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    cp = get_default_config()
    cs = get_configstr(cp)

    logging.debug(f"got config: {cs}")

    if args.setup:
        s = SetUp(cp)
        s.execute()

    if args.assign:
        dq = Queue()
        for solooutdir in args.soloutdirs:
            fq = AssignCellType(cp, solooutdir, outlist)
            dq.put(fq)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('sra', 'max_downloads'))
        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')
