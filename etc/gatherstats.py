# aggregates statistics.

import pandas as pd
from scipy.io import mmread
import os
import scanpy as sc  # pip install
import numpy as np
import io
import argparse

import subprocess
import sys

from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty

from scqc.utils import *

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


class GetStats(object):

    def __init__(self, config, srpid, ):
        self.log = logging.getLogger('stats')
        self.config = config
        self.solooutdir = srpid    # Solo.out/
        self.staroutdir = os.path.expanduser(
            self.config.get('stats', 'staroutdir'))
        self.statdir = os.path.expanduser(self.config.get('stats', 'statdir'))

        self.metadir = os.path.expanduser(self.config.get('stats', 'metadir'))

    def _gather_stats_from_STAR(self):

        with open(f"{self.solooutdir}/Barcodes.stats", 'r') as f:
            lines = f.readlines()
            lines = [" ".join(line.split()) for line in lines]
            lines = [line.split() for line in lines]
            barcode_stats = pd.DataFrame(lines, columns=["stat", "value"])

        with open(f"{self.solooutdir}/Gene/Features.stats", 'r') as f:
            lines = f.readlines()
            lines = [" ".join(line.split()) for line in lines]
            lines = [line.split() for line in lines]
            feature_stats = pd.DataFrame(lines, columns=["stat", "value"])

        summary_stats = pd.read_csv(
            f"{self.solooutdir}/Gene/Summary.csv", sep=",", header=None)
        summary_stats.columns = ["stat", "value"]

        acc = self.solooutdir.split("/")[-1].split("Solo.out")[0][:-1]
        barcode_stats['stat_source'] = 'barcode'
        feature_stats['stat_source'] = 'feature'
        summary_stats['stat_source'] = 'summary'
        barcode_stats["accession"] = acc
        feature_stats["accession"] = acc
        summary_stats["accession"] = acc

        return(barcode_stats, feature_stats, summary_stats)

    def _parse_STAR_mtx(self):
        # note that scanpy uses cell x gene.
        # read into anndata
        # path should end with "Solo.out"
        # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"
        path = f'{self.solooutdir}/Gene/filtered'
        # mtx_files = os.listdir(path)
        adata = sc.read(f'{path}/matrix.mtx').T

        genenames = pd.read_csv(f'{path}/features.tsv',
                                sep="\t", header=None, dtype="str")
        genenames.columns = ['gene_accession', 'gene_symbol', 'source']
        # genenames.index = genenames.gene_accession

        cellids = pd.read_csv(f'{path}/barcodes.tsv', sep="\t", header=None)
        cellids.columns = ["cell_id"]
        # cellids.index = cellids.cell_id

        adata.var = genenames
        adata.obs = cellids

        return adata

    def _get_stats_scanpy(self, adata):
        adata.var['mt'] = adata.var.gene_symbol.str.startswith('mt-')
        adata.var['ERCC'] = adata.var.gene_symbol.str.startswith('ERCC')
        # adata should be in raw counts
        sc.pp.calculate_qc_metrics(
            adata, expr_type='counts', var_type='genes', qc_vars=['mt', "ERCC"],
            percent_top=(50, 100, 200, 500), inplace=True, use_raw=False)

        return adata

    def execute(self):
        # outdir = "/home/johlee/scqc/stats"
        # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"
        barcode_stats, feature_stats, summary_stats = self._gather_stats_from_STAR(
            self.solooutdir)

        all_stats = summary_stats.append(barcode_stats).append(feature_stats)
        # did we already save these?
        acc = barcode_stats.accession.unique()[0]
        fname = f'{self.outdir}/{acc}_starout_stats.tsv'
        runs_done = pd.read_csv(
            f'{self.solooutdir}/Gene/raw/barcodes.tsv', sep="\t", header=None)
        projectout = f'{self.outdir}/{acc}_runs_done.tsv'

        all_stats.to_csv(fname, sep="\t", index=None,
                         header=not os.path.isfile(fname), mode='w')
        runs_done.to_csv(projectout, sep="\t", index=None,
                         header=False, mode='a')

        adata = self._parse_STAR_mtx(self.solooutdir)
        adata = self._get_stats_scanpy(adata)

        # scanpy default is to overwrite existing file
        h5file = f'{self.outdir}/{acc}.h5ad'
        adata.write(h5file)


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

    parser.add_argument('-s', '--srp',
                        metavar='srpid',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='')



    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    cp = get_default_config()
    cs = get_configstr(cp)

    logging.debug(f"got config: {cs}")

    if args.fasterq is not None:
        dq = Queue()
        for srr in args.fasterq:
            fq = FasterqDump(cp, srr)
            dq.put(fq)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('sra', 'max_downloads'))
        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')
