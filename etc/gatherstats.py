# aggregates statistics.

import argparse
import io
import logging
import os
import subprocess
import sys
from configparser import ConfigParser
from queue import Empty, Queue
from threading import Thread

import numpy as np
import pandas as pd
import scanpy as sc  # pip install
from scipy.io import mmread
from scqc.utils import gini_coefficient, sparse_pairwise_corr

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
        self.starindexdir = os.path.expanduser(
            self.config.get('stats', 'starindexdir'))
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

    # need to pass in star index directory

    def _parse_STAR_mtx(self):
        # note that scanpy uses cell x gene.
        # read into anndata
        # path should end with "Solo.out"
        # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"
        path = f'{self.solooutdir}/Gene/filtered'

        # mtx_files = os.listdir(path)
        adata = sc.read(f'{path}/matrix.mtx').T

        geneinfo = pd.read_csv(
            '~/scqc/supplement_data/genomes/mouse/STAR_index/geneInfo.tab', sep="\t",  skiprows=1, header=None, dtype=str)
        geneinfo.columns = ['gene_accession', 'gene_symbol', 'type']
        # geneinfo.index = geneinfo.gene_accession

        genenames = pd.read_csv(f'{path}/features.tsv',
                                sep="\t", header=None, dtype=str)
        genenames.columns = ['gene_accession', 'gene_symbol', 'source']

        genenames = genenames.merge(geneinfo, how='left', on=[
            "gene_accession", 'gene_symbol'], indicator=True)
        genenames.index = genenames.gene_accession

        cellids = pd.read_csv(f'{path}/barcodes.tsv', sep="\t", header=None)
        cellids.columns = ["cell_id"]
        # cellids.index = cellids.cell_id

        adata.var = genenames
        adata.obs = cellids

        return adata

    def _get_stats_scanpy(self, adata):

        # consider different gene sets - ERCC  corresponds to spike ins
        adata.var['mt'] = adata.var.gene_symbol.str.startswith('mt-')
        adata.var['ERCC'] = adata.var.gene_symbol.str.startswith('ERCC')
        adata.var['ribo'] = adata.var.type == "rRNA"
        adata.var['cytoplasm'] = None       # GO Term/kegg?
        adata.var['metabolism'] = None      # GO Term/kegg?
        adata.var['membrane'] = None        # GO Term/kegg?

        qcvars = ['mt', 'ERCC', 'ribo', 'female', 'male',
                  'essential', 'cell_cycle']

        # get housekeeping genes
        hk_genes = pd.read_csv(self.housekeeping, sep=",")
        adata.var['housekeeping'] = adata.var.gene_symbol.isin(
            hk_genes.gene)

        # get gender markers
        female_genes = pd.read_csv(self.female_markers, sep=",")
        adata.var['female'] = adata.var.gene_symbol.isin(
            female_genes.gene)

        male_genes = pd.read_csv(self.male_markers, sep=",")
        adata.var['male'] = adata.var.gene_symbol.isin(
            male_genes.gene)

        # get essential genes
        essential = pd.read_csv(self.essential, sep=",")
        adata.var['essential'] = adata.var.gene_symbol.isin(
            essential.gene)

        # get cell cycle genes
        cc = pd.read_csv(self.cc_marker_path, sep=",")
        adata.var['cell_cycle'] = adata.var.gene_symbol.isin(cc.gene)
        for i in cc.cluster.unique():
            adata.var[f'cc_cluster_{i}'] = adata.var.cluster == i
            qcvars.append(f'cc_cluster_{i}')

        # append stats for gene sets above using scanpy
        sc.pp.calculate_qc_metrics(
            adata,
            expr_type='counts', var_type='genes',
            percent_top=(50, 100, 200, 500), inplace=True, use_raw=False,
            qc_vars=qcvars)

        # computes the N+M x N+M corrcoef matrix - extract off diagonal block
        adata.obs['corr_to_mean'] = np.array(sparse_pairwise_corr(
            adata.var.mean_counts, adata.X)[0, 1:]).flatten()

        # unstructured data - dataset specific
        adata.uns['gini_by_counts'] = gini_coefficient(
            adata.obs['total_counts'])

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
    pass
