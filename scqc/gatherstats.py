#!/usr/bin/env python


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
from sklearn.metrics.pairwise import euclidean_distances
# import matplotlib.cbook as cbook


# SRP090110 - mouse brain SS (500 cells)    Not found in scbrain_pids_mouse!
# SRP124513 - mouse brain SS (1700 cells)   Not found in scbrain_pids_mouse!
# SRP110034 - mouse brain 10x (1700 cells)  Not found in scbrain_pids_mouse!
# SRP106908 - mouse brain SS (3186 cells) (atlas)
# SRP135960 - mouse brain 10x (509000 cells) (atlas)

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *

sc.settings.verbosity = 4


LOGLEVELS = {
    10: 'debug',
    20: 'info',
    30: 'warn',
    40: 'err',
    50: 'fatal',
}

class GetStats(object):

    def __init__(self, config):
        self.log = logging.getLogger('stats')
        self.config = config
        self.species = self.config.get('stats', 'species')
        self.outputdir = os.path.expanduser(
            self.config.get('stats', 'outputdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('stats', 'resourcedir'))
        self.metadir = os.path.expanduser(
            self.config.get('stats', 'metadir'))
        self.tempdir = os.path.expanduser(
            self.config.get('stats', 'tempdir'))
        self.cachedir = os.path.expanduser(
            self.config.get('stats', 'cachedir'))
        # gene sets

        self.cell_cycle_genes = os.path.expanduser(
            self.config.get('stats', 'cell_cycle_genes'))
        self.stable_housekeepinig = os.path.expanduser(
            self.config.get('stats', 'stable_housekeepinig'))
        self.essential_genes = os.path.expanduser(
            self.config.get('stats', 'essential_genes'))
        self.female_genes = os.path.expanduser(
            self.config.get('stats', 'female_genes'))
        self.male_genes = os.path.expanduser(
            self.config.get('stats', 'male_genes'))


        # hvg params 
        self.hvg_min_mean = self.config.get('stats', 'hvg_min_mean')
        self.hvg_max_mean = self.config.get('stats', 'hvg_max_mean')
        self.hvg_min_disp = self.config.get('stats', 'hvg_min_disp')
        self.hvg_max_disp = self.config.get('stats', 'hvg_max_disp')
        self.hvg_flavor = self.config.get('stats', 'hvg_flavor')
        #pca params
        self.nPCA = self.config.get('stats', 'nPCA')
        self.use_hvg = self.config.get('stats', 'use_hvg')
        # neighbor graph params
        self.n_neighbors = self.config.get('stats', 'n_neighbors')

        # metamarker params
        self.class_markerset =  os.path.expanduser(
            self.config.get('metamarker', 'class_markerset'))
        self.subclass_markerset =  os.path.expanduser(
            self.config.get('metamarker', 'subclass_markerset'))
        self.max_rank = self.config.get('metamarker', 'max_rank')


    def execute(self, srpid):
        '''
        STAR outputs should be placed in the temp directory 
        Look through all of the solo out directories for the project,
        aggregate them, and run stats collectively

        '''

        # outdir = "/home/johlee/scqc/stats"
        solooutdirs = f'{self.cachedir}/{srpid}'
        solooutdirs = glob.glob(f'{solooutdirs}/*Solo.out')

        self.log.debug(f'starting with projectid {srpid}')

        adatas = []   # loops through all star runs in the project
        for solooutdir in solooutdirs:
            self.log.debug(f'opening {solooutdir}')
            # open one soloutdir, store as adata. 
            adatas.append(self._parse_STAR_mtx(solooutdir))

        # merge all of the datasets 
        # technology is obtained from the soloout dir name
        batchnames = [ os.path.basename(solooutdir).replace('_Solo.out','') for i in range(len(solooutdirs)) ]
        ids = pd.DataFrame([ batchname.split('_') for batchname in batchnames] )
        self.log.debug(f'starting with projectid {ids}')

        ids.columns = ['id','tech']
        # which are runs and which are projects? 
        ids['run_id'] = ids.id[ids.id.str.contains('RR')]
        ids['proj_id'] = ids.id[ids.id.str.contains('RP')]
        
        # fill in proj_id NaNs
        ids.proj_id[ids.proj_id.isnull()] = srpid

        self.log.debug(f'joining adatas for {srpid}')
        adata = adatas[0].concatenate(adatas[1:], batch_categories = ids.id )
        adata.obs.columns = ['cell_id','id']
        ind = adata.obs.index

        tmpdf = pd.merge(adata.obs, ids, on= 'id')   # include technology in obs.
        tmpdf.index = tmpdf.cell_id
        adata.obs = tmpdf 
        adata.obs.index= ind

        # fill in missing run_ids with the cell_id. (typically smart seq experiments)
        adata.obs.run_id[adata.obs.run_id.isnull()] = adata.obs.cell_id[adata.obs.run_id.isnull()]

        # get batch information
        self.log.debug(f'getting batch info for {srpid}')

        impute = pd.read_csv(f'{self.metadir}/impute.tsv', sep='\t',index_col=0)
        impute = impute.loc[impute.proj_id == srpid,: ]

        tmpdf =  pd.merge(adata.obs, impute )[['cell_id','run_id','exp_id', 'samp_id','proj_id','tech', 'batch' ]]
        tmpdf.index= ind

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

        self.log.debug(f'assinging cell types using metamarkers for {srpid}')
        tmp_path = self._run_MetaMarkers(h5file)
        self.log.debug(f'done with metamarkers for {srpid} - saving to h5file')
        
        tmpdf = pd.read_csv(tmp_path ,sep="\t", index_col =0)
        self.log.debug(f'read df')


        # TODO separate into the scores - enrichment - prediction
        tmp =  adata.obs.merge(tmpdf, left_on = 'cell_id' ,right_on = 'cell')
        tmp.index = tmp.cell_id
        tmp.index.name = None
        tmp.reindex(adata.obs.index)
        adata.obsm['MetaMarkers'] = tmp

        self.log.debug('Saving to h5ad')
        adata.write(h5file)
    
    # currently not used 
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
        adata = sc.read_mtx(f'{path}/matrix.mtx').T

        
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
        # adata.obs['batch'] = None

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
        # adata.uns['gini_by_counts'] = gini_coefficient(
        #     adata.obs['total_counts'])

        self.log.debug('computing hvg')
        # highly variable genes - expects log data
        sc.pp.highly_variable_genes(adata,
                                    min_mean=float(self.hvg_min_mean),
                                    max_mean=float(self.hvg_max_mean),
                                    min_disp=float(self.hvg_min_disp),
                                    max_disp=float(self.hvg_max_disp),
                                    flavor=self.hvg_flavor)

        # umap coords

        ncomp = min(adata.shape[0] -1 , int(self.nPCA) )
        self.log.debug('running pca\n')
        sc.tl.pca(adata,
                  n_comps=ncomp, zero_center=False, use_highly_variable=self.use_hvg)
        print('getting neighbors\n')
        sc.pp.neighbors(adata, n_neighbors=int(self.n_neighbors), n_pcs=ncomp)
        print('getting umap\n')
        sc.tl.umap(adata)       # umap coords are in adata.obsm['X_umap']
        
        adata.obs.cell_id = adata.obs.index
        return adata

    def _run_EGAD_by_batch(self, adata, rank_standardized = False, 
        batch_column = ['run_id','exp_id','samp_id','proj_id','batch','class_predicted','subclass_predicted'] ):

        # XXX memory intensive for datasets with large number of cells.
        # build the pairwise distance matrix

        nw = euclidean_distances(adata.X)
        if rank_standardized:
            nw = rank(nw.max()-nw )
        else :
            nw = 1- nw / nw.max()

        nw = pd.DataFrame( nw )
    
        
        # convert to a cell x batch binary df
        go = pd.DataFrame( )
        
        # for exp_id in adata.obs.exp_id.unique() :
        #     go[exp_id] = adata.obs.exp_id == exp_id
        for col in batch_column:
            for id in adata.obs[col].unique() :
                go[id] = adata.obs[col] == id
        go = go.reset_index(drop=True)   
        # batches should contain as least 10 cells
        # and no more that 90%  of all cells 
        go = go.T.drop_duplicates().T
        res = run_egad(go, nw,  nFold=3, min_count=5, max_count=np.ceil(adata.shape[0] * .95 ) )
        
        adata.uns['EGAD'] = res
        return( adata)
    
    def _run_MetaMarkers(self, h5path):
        '''
        h5path should be saved in the temp directory before this stage. 
        After MetaMarkers, move completed h5ad file to output/
        '''

        scriptpath = '/home/johlee/git/scqc/scqc/get_markers.R'    # should already be in path
        # should have only one outpath per dataset. 
        # TODO verify file does not already exist. Think about what to do if it does
        # outpath = f'{self.outputdir}/{os.path.basename(h5path)}'
        cmd = ['Rscript', f'{scriptpath}',
               '--class_markerset', f'{self.class_markerset}',
               '--subclass_markerset', f'{self.subclass_markerset}',
               '--h5path', f'{h5path}',
               '--max_rank', f'{self.max_rank}'
               ] #'--outprefix', f'{outpath}'
        subprocess.run(cmd)

        tmp_path = h5path.replace('.h5ad','.tsv')
        return(tmp_path)
    

# class MetaMarker(object):
    
#     def __init__(self, config):
#         self.log = logging.getLogger('statistics')
#         self.config = config
#         self.species = self.config.get('statistics', 'species')
#         self.outputdir = os.path.expanduser(
#             self.config.get('statistics', 'outputdir'))
#         self.resourcedir = os.path.expanduser(
#             self.config.get('statistics', 'resourcedir'))
#         self.metadir = os.path.expanduser(
#             self.config.get('statistics', 'metadir'))
#         self.tempdir = os.path.expanduser(
#             self.config.get('statistics', 'tempdir'))

#         # metamarker params
#         self.class_markerset =  os.path.expanduser(
#             self.config.get('metamarker', 'class_markerset'))
#         self.subclass_markerset =  os.path.expanduser(
#             self.config.get('metamarker', 'subclass_markerset'))
#         self.max_rank = self.config.get('metamarker', 'max_rank')


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
