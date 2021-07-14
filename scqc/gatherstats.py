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
import matplotlib.pyplot as plt
import seaborn as sns       # inputs as dataframes
import sklearn
# import matplotlib.cbook as cbook


# SRP090110 - mouse brain SS (500 cells)    Not found in scbrain_pids_mouse!
# SRP124513 - mouse brain SS (1700 cells)   Not found in scbrain_pids_mouse!
# SRP110034 - mouse brain 10x (1700 cells)  Not found in scbrain_pids_mouse!
# SRP106908 - mouse brain SS (3186 cells) (atlas)
# SRP135960 - mouse brain 10x (509000 cells) (atlas)

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

class GetStats(object):

    def __init__(self, config):
        self.log = logging.getLogger('statistics')
        self.config = config
        self.species = self.config.get('statistics', 'species')
        self.outputdir = os.path.expanduser(
            self.config.get('statistics', 'outputdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('statistics', 'resourcedir'))
        self.metadir = os.path.expanduser(
            self.config.get('statistics', 'metadir'))
        self.tempdir = os.path.expanduser(
            self.config.get('statistics', 'tempdir'))

        # gene sets

        self.cell_cycle_genes = os.path.expanduser(
            self.config.get('statistics', 'cell_cycle_genes'))
        self.stable_housekeepinig = os.path.expanduser(
            self.config.get('statistics', 'stable_housekeepinig'))
        self.essential_genes = os.path.expanduser(
            self.config.get('statistics', 'essential_genes'))
        self.female_genes = os.path.expanduser(
            self.config.get('statistics', 'female_genes'))
        self.male_genes = os.path.expanduser(
            self.config.get('statistics', 'male_genes'))


        # hvg params 
        self.hvg_min_mean = self.config.get('statistics', 'hvg_min_mean')
        self.hvg_max_mean = self.config.get('statistics', 'hvg_max_mean')
        self.hvg_min_disp = self.config.get('statistics', 'hvg_min_disp')
        self.hvg_max_disp = self.config.get('statistics', 'hvg_max_disp')
        self.hvg_flavor = self.config.get('statistics', 'hvg_flavor')
        #pca params
        self.nPCA = self.config.get('statistics', 'nPCA')
        self.use_hvg = self.config.get('statistics', 'use_hvg')
        # neighbor graph params
        self.n_neighbors = self.config.get('statistics', 'n_neighbors')
        self.n_pcs = self.config.get('statistics', 'n_pcs')



    def execute(self, srpid):
        # outdir = "/home/johlee/scqc/stats"
        solooutdirs = f'{self.outputdir}/{srpid}'
        solooutdirs = glob.glob(f'{solooutdirs}/*Solo.out')

        self.log.debug(f'starting with projectid {srpid}')

        adatas = []   # loops through all star runs in the project
        for solooutdir in solooutdirs:
            self.log.debug(f'opening {solooutdir}')
            # open one soloutdir, store as adata. 
            adatas.append(self._parse_STAR_mtx(solooutdir))

        # merge all of the datasets 
        batchnames = [ os.path.basename(solooutdirs[i]).replace('_Solo.out','') for i in range(len(solooutdirs)) ]
        ids = pd.DataFrame([ batchname.split('_') for batchname in batchnames] )

        ids.columns = ['id','tech']
        # which are runs and which are projects? \
        ids['run_id'] = ids.id[ids.id.str.contains('RR')]
        ids['proj_id'] = ids.id[ids.id.str.contains('RP')]
        
        # fill in proj_id NaNs
        ids.proj_id[ids.proj_id.isnull()] = srpid

        self.log.debug(f'joining adatas for {srpid}')
        adata = adatas[0].concatenate(adatas[1:], batch_categories = ids.id )
        adata.obs.batch='id'
        
        tmpdf = pd.merge(adata.obs, ids, on= 'id')   # include technology in obs.
        tmpdf.index = tmpdf.cell_id
        tmpdf.index.name =None
        adata.obs = tmpdf 

        # fill in missing run_ids with the cell_id. (typically smart seq experiments)
        adata.obs.run_id[adata.obs.run_id.isnull()] = adata.obs.cell_id[adata.obs.run_id.isnull()]

        # get batch information
        self.log.debug(f'getting batch info for {srpid}')

        impute = pd.read_csv(f'{self.metadir}/impute.tsv', sep='\t',index_col=0)
        impute = impute.loc[impute.proj_id == srpid,: ]

        tmpdf =  pd.merge(adata.obs, impute )[['cell_id','run_id','exp_id', 'samp_id','proj_id','tech', 'batch' ]]
        tmpdf.index = tmpdf.cell_id
        tmpdf.index.name =None
        adata.obs = tmpdf

        adata = self._get_stats_scanpy(adata)

        # XXX scanpy default is to overwrite existing file
        # basename = os.path.basename(solooutdir)
        # basename = basename.replace('_Solo.out', '.h5ad')
        h5file = f'{self.outputdir}/{srpid}.h5ad'
        adata.var.index.name = None
        self.log.debug(adata.obs)
        self.log.debug(adata.var)

        # for testing...
        # adata[solooutdir].obs.to_csv(f'{self.tempdir}/obs.tsv',sep="\t")
        # adata[solooutdir].var.to_csv(f'{self.tempdir}/var.tsv',sep="\t")
        # consolidate these...
        self.log.debug(f'saving to {self.outputdir}/{srpid}.h5ad ')

        adata.write(h5file)

    def _gather_stats_from_STAR(self, solooutdir):

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

    def _parse_STAR_mtx(self, solooutdir):
        # note that scanpy uses cell x gene.
        # read into anndata
        # path should end with "Solo.out"
        # solooutdir = "/home/johlee/scqc/starout/SRP308826_smartseq_Solo.out"
        path = f'{solooutdir}/Gene/filtered'

        # mtx_files = os.listdir(path)
        adata = sc.read(f'{path}/matrix.mtx').T

        
        geneinfo = pd.read_csv(f'{self.resourcedir}/{self.species}/geneInfo.tab', sep="\t",  skiprows=1, header=None, dtype=str)
        geneinfo.columns = ['gene_accession', 'gene_symbol', 'type']
        # geneinfo.index = geneinfo.gene_accession
        # geneinfo.index.name =None
        genenames = pd.read_csv(f'{path}/features.tsv',
                                sep="\t", header=None, dtype=str)
        genenames.columns = ['gene_accession', 'gene_symbol', 'source']

        genenames = genenames.merge(geneinfo, how='left', on=[
            "gene_accession", 'gene_symbol'], indicator=True)
        genenames.index = genenames.gene_accession
        genenames.index.name =None
        cellids = pd.read_csv(f'{path}/barcodes.tsv', sep="\t", header=None)
        cellids.columns = ["cell_id"]

        cellids.index = cellids.cell_id
        cellids.index.name =None

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
        # adata.var['cytoplasm'] = None       # GO Term/kegg?
        # adata.var['metabolism'] = None      # GO Term/kegg?
        # adata.var['membrane'] = None        # GO Term/kegg?

        qcvars = ['mt', 'ribo', 'female', 'male',
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
        self.log.debug(f"ess:{essential}")
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
        
        # natural log thge data
        sc.pp.log1p(adata)  

        # batch this.
        # X = np.random.normal(size=(10, 3))
        # F = np.zeros((10, ))

        # import multiprocessing
        # pool = multiprocessing.Pool(processes=16)
        # if number of processes is not specified, it uses the number of core
        # F[:] = pool.map(my_function, (X[i,:] for i in range(10)) )
        tmp = gini_coefficient_fast(adata.X)

        adata.obs['gini']  = tmp
        # unstructured data - dataset specific
        adata.uns['gini_by_counts'] = gini_coefficient(
            adata.obs['total_counts'])

        self.log.debug('computing hvg')
        # highly variable genes - expects log data
        sc.pp.highly_variable_genes(adata,
                                    min_mean=float(self.hvg_min_mean),
                                    max_mean=float(self.hvg_max_mean),
                                    min_disp=float(self.hvg_min_disp),
                                    max_disp=float(self.hvg_max_disp),
                                    flavor=self.hvg_flavor)

        # umap coords
        sc.tl.pca(adata,
                  n_comps=int(self.nPCA), zero_center=False, use_highly_variable=self.use_hvg)

        sc.pp.neighbors(adata, n_neighbors=int(self.n_neighbors), n_pcs=int(self.n_pcs))
        sc.tl.umap(adata)       # umap coords are in adata.obsm['X_umap']

        return adata

    def run_EGAD_by_batch(self,adata):

        # XXX memory intensive for datasets with large number of cells.
        # build the pairwise distance matrix
        nw = sklearn.metrics.pairwise.euclidean_distances(adata.X)
        nw = pd.DataFrame(1- nw / nw.max() )
    
        
        # convert to a cell x batch binary df
        go = pd.DataFrame( )
        # for exp_id in adata.obs.exp_id.unique() :
        #     go[exp_id] = adata.obs.exp_id == exp_id
        for samp_id in adata.obs.samp_id.unique() :
            go[samp_id] = adata.obs.samp_id == samp_id
        for run_id in adata.obs.run_id.unique() :
            go[run_id] = adata.obs.run_id == run_id
        # batches should contain as least 10 cells
        # and no more that 75%  of all cells 
        res = run_egad(go, nw,  nFold=3, min_count=10, max_count=np.ceil(adata.shape[0] * .75 ) )
        

if __name__ == "__main__":

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
                        help='gather stats')

    parser.add_argument('-f', '--figures',
                        metavar='projectids',
                        type=str,
                        nargs='+',
                        default=None,
                        help='build figures')

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
