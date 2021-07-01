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
import glob 
import numpy as np
import pandas as pd
import scanpy as sc  # pip install
from scipy.io import mmread
from scipy.sparse import base


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *

LOGLEVELS = {
    10: 'debug',
    20: 'info',
    30: 'warn',
    40: 'err',
    50: 'fatal',
}

TO_PLOT = ['reads']

#TODO given just proj_id, split to smart seq and 10x and proceed independently.
#TODO generate figures for the dataset
class GetStats(object):

    def __init__(self, config ):
        self.log = logging.getLogger('statistics')
        self.config = config
        self.outputdir = os.path.expanduser(
            self.config.get('statistics', 'outputdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('statistics', 'resourcedir'))
        self.metadir = os.path.expanduser(self.config.get('statistics', 'metadir'))
        self.tempdir = os.path.expanduser(self.config.get('statistics', 'tempdir'))

        # gene sets
        self.cell_cycle_genes = os.path.expanduser(self.config.get('statistics', 'cell_cycle_genes'))
        self.stable_housekeepinig = os.path.expanduser(self.config.get('statistics', 'stable_housekeepinig'))
        self.essential_genes = os.path.expanduser(self.config.get('statistics', 'essential_genes'))
        self.female_genes = os.path.expanduser(self.config.get('statistics', 'female_genes'))
        self.male_genes = os.path.expanduser(self.config.get('statistics', 'male_genes'))

        self.nPCA = self.config.get('statistics', 'nPCA')

    def execute(self,srpid):
        # outdir = "/home/johlee/scqc/stats"
        solooutdirs = f'{self.outputdir}/{srpid}'
        solooutdirs = glob.glob(f'{solooutdirs}/*Solo.out')

        self.log.debug(f'starting with projectid {srpid}')

        adata ={}
        for solooutdir in solooutdirs :
            self.log.debug(f'looking at {solooutdir}')
            # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"

            # TODO fix this!  
            # Dumps stats from STAR into one file and stores done list

            # barcode_stats, feature_stats, summary_stats = self._gather_stats_from_STAR(
            #     solooutdir)
            # all_stats = summary_stats.append(barcode_stats).append(feature_stats)
            # # did we already save these?
            # acc = barcode_stats.accession.unique()[0]
            # fname = f'{self.metadata}/starout_stats.tsv'
            # self.log.debug(f'read {fname}')

            # runs_done = pd.read_csv(
            #     f'{solooutdir}/Gene/raw/barcodes.tsv', sep="\t", header=None)
            # projectout = f'{self.outputdir}/{acc}_runs_done.tsv'

            # all_stats.to_csv(fname, sep="\t", index=None,
            #                 header=not os.path.isfile(fname), mode='w')
            # runs_done.to_csv(projectout, sep="\t", index=None,
            #                 header=False, mode='a')

            adata[solooutdir] = self._parse_STAR_mtx(solooutdir)
            adata[solooutdir] = self._get_stats_scanpy(adata[solooutdir])

            # scanpy default is to overwrite existing file
            basename = os.path.basename(solooutdir)
            basename =basename.replace('_Solo.out' ,'.h5ad')
            h5file = f'{self.outputdir}/{srpid}/{basename}'
            adata[solooutdir].var.index.name = None
            self.log.debug(adata[solooutdir].obs.columns)
            self.log.debug(adata[solooutdir].var.columns)
            

            # for testing...
            adata[solooutdir].obs.to_csv(f'{self.tempdir}/obs.tsv',sep="\t")
            adata[solooutdir].var.to_csv(f'{self.tempdir}/var.tsv',sep="\t")
            # consolidate these... 
            adata[solooutdir].write(h5file)



    def _gather_stats_from_STAR(self,solooutdir):

        with open(f"{solooutdir}/Barcodes.stats", 'r') as f:
            lines = f.readlines()
            lines = [" ".join(line.split()) for line in lines]
            lines = [line.split() for line in lines]
            barcode_stats = pd.DataFrame(lines, columns=["stat", "value"])

        with open(f"{solooutdir}/Gene/Features.stats", 'r') as f:
            lines = f.readlines()
            lines = [" ".join(line.split()) for line in lines]
            lines = [line.split() for line in lines]
            feature_stats = pd.DataFrame(lines, columns=["stat", "value"])

        summary_stats = pd.read_csv(
            f"{solooutdir}/Gene/Summary.csv", sep=",", header=None)
        summary_stats.columns = ["stat", "value"]

        acc = solooutdir.split("/")[-1].split("Solo.out")[0][:-1]
        barcode_stats['stat_source'] = 'barcode'
        feature_stats['stat_source'] = 'feature'
        summary_stats['stat_source'] = 'summary'
        barcode_stats["accession"] = acc
        feature_stats["accession"] = acc
        summary_stats["accession"] = acc

        return(barcode_stats, feature_stats, summary_stats)

    # need to pass in star index directory

    def _parse_STAR_mtx(self,solooutdir):
        # note that scanpy uses cell x gene.
        # read into anndata
        # path should end with "Solo.out"
        # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"
        path = f'{solooutdir}/Gene/filtered'

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

        # gather imputed data and append to obs

        return adata

    def _get_stats_scanpy(self, adata):
        adata.obs['batch'] = None


        # consider different gene sets - ERCC  corresponds to spike ins
        adata.var['mt'] = adata.var.gene_symbol.str.startswith('mt-')
        # adata.var['ERCC'] = adata.var.gene_symbol.str.startswith('ERCC')
        adata.var['ribo'] = adata.var.type == "rRNA"
        adata.var['cytoplasm'] = None       # GO Term/kegg?
        adata.var['metabolism'] = None      # GO Term/kegg?
        adata.var['membrane'] = None        # GO Term/kegg?

        qcvars = ['mt','ribo', 'female', 'male',
                  'essential', 'cell_cycle']

        # get housekeeping genes
        self.log.debug(f"gathering marker sets:")
        hk_genes = pd.read_csv(self.stable_housekeepinig, sep=",")
        self.log.debug(f"housekeeping:{hk_genes}")
        adata.var['housekeeping'] = adata.var.gene_symbol.isin(
            hk_genes.gene)


        # get gender markers
        female_genes = pd.read_csv(self.female_genes, sep=",")
        self.log.debug(f"female:{female_genes}")
        adata.var['female'] = adata.var.gene_symbol.isin(
            female_genes.gene)

        male_genes = pd.read_csv(self.male_genes, sep=",")
        self.log.debug(f"male:{male_genes}")
        adata.var['male'] = adata.var.gene_symbol.isin(
            male_genes.gene)

        # get essential genes
        essential = pd.read_csv(self.essential_genes, sep=",")
        self.log.debug(f"ess:{essential}" )
        adata.var['essential'] = adata.var.gene_symbol.isin(
            essential.gene)

        # get cell cycle genes
        cc = pd.read_csv(self.cell_cycle_genes, sep=",")
        self.log.debug(f"cc:{cc}")
        adata.var['cell_cycle'] = adata.var.gene_symbol.isin(cc.gene)
        for i in cc.cluster.unique():
            adata.var[f'cc_cluster_{i}'] = adata.var.cell_cycle == i
            qcvars.append(f'cc_cluster_{i}')

        self.log.debug('calculating qc metrics')
        sc.pp.calculate_qc_metrics(
            adata,
            expr_type='counts', var_type='genes',
            percent_top=(50, 100, 200, 500), inplace=True, use_raw=False,
            qc_vars=qcvars)

        # computes the N+M x N+M corrcoef matrix - extract off diagonal block
        self.log.debug('calculating corrtomean')
        adata.obs['corr_to_mean'] = np.array(sparse_pairwise_corr(
            adata.var.mean_counts, adata.X)[0, 1:]).flatten()

        # this part is slow (~10-15 min) because of a for loop - can parallelize
        self.log.debug('calculating gini for each cell')
        # NOTE SLOW commented for testing
        adata.obs['gini'] = gini_coefficient_spmat(adata.X )

        # unstructured data - dataset specific
        adata.uns['gini_by_counts'] = gini_coefficient(
            adata.obs['total_counts'])

        # highly variable genes - expects log data
        sc.pp.log1p(adata) # natural log
        sc.pp.highly_variable_genes( adata,
                min_mean=self.hvg_min_mean,
                max_mean=self.hvg_max_mean,
                min_disp=self.hvg_min_disp,
                max_disp=self.hvg_max_disp)


        # XXX  does it make sense to have a umap for every run? Or one for the project? 
        # we'll do both for now
        # umap coords
        sc.tl.pca(adata,
             n_comps = self.nPCA, zero_center=False,use_highly_variable=self.use_hvg)
        
        sc.pp.neighbors(adata, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
        sc.tl.umap(adata)       # umap coords are in adata.obsm['X_umap']
        

        return adata

    # consolidate h5ads for the project.
    

class BuildFigures(object):
    def __init__():
        pass

    

    pass

if __name__ =="__main__":

    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    # logging.getLogger().setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug',
                        action="store_true",
                        dest='debug',
                        help='debug logging')

    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        dest='verbose',
                        help='verbose logging')

            
    parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default='~/git/scqc/etc/scqc.conf',
                            help='Config file path [~/git/scqc/etc/scqc.conf]')

    parser.add_argument('-p', '--projectids',
                        metavar='projectids',
                        type=str,
                        nargs='+',
                        default=None,
                        help='Perform tech and batch imputation on supplied projectid')


    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        cp = get_default_config()
       
    cs = get_configstr(cp)
    logging.debug(f"got config: {cs}")

    logging.debug(f"args: {args}")

    if args.projectids is not None:
        q = GetStats(cp)
        print(args.projectids)
        for pid in args.projectids:
        # args.projectid is a list of project ids
            q.execute(pid)
