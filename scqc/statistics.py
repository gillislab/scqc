#!/usr/bin/env python
#
# Aggregates statistics. Outputs to /output/<proj_id>.h5ad
#

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
import scanpy as sc  
from scipy.io import mmread
from scipy.sparse import base
import matplotlib.pyplot as plt
import seaborn as sns       # inputs as dataframes
from sklearn.metrics.pairwise import euclidean_distances
# import matplotlib.cbook as cbook



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

class Statistics(object):

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
        self.cachedir = os.path.expanduser(
            self.config.get('statistics', 'cachedir'))
        self.nocleanup = self.config.getboolean('statistics','nocleanup')
        
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
        
        # metamarker params
        self.class_markerset =  os.path.expanduser(
            self.config.get('metamarker', 'class_markerset'))
        self.subclass_markerset =  os.path.expanduser(
            self.config.get('metamarker', 'subclass_markerset'))
        self.max_rank = self.config.get('metamarker', 'max_rank')
        



    def execute(self, proj_id):
        '''
        STAR outputs should be placed in the temp directory 
        Look through all of the solo out directories for the project,
        aggregate them, and run stats collectively

        '''
        self.log.debug(f'starting project id {proj_id}')
        done = None
        part = None
        seen = proj_id

        try:
            solooutdirs = glob.glob(f'{self.cachedir}/{proj_id}/*Solo.out')
            adatas = []   # loops through all star runs in the project
            for solooutdir in solooutdirs:
                self.log.debug(f'opening {solooutdir}')
                # open one soloutdir, store as adata. 
                adatas.append(self._parse_STAR_mtx(solooutdir))
    
            # merge all of the datasets 
            # technology is obtained from the soloout dir name
            batchnames = [ os.path.basename(dir).replace('_Solo.out','') for dir in solooutdirs ]
            ids = pd.DataFrame([ batchname.split('_') for batchname in batchnames] )
            self.log.debug(f'starting with projectid {ids}')
    
            ids.columns = ['id','tech']
            # which are runs and which are projects? 
            ids['run_id'] = ids.id[ids.id.str.contains('RR')]
            ids['proj_id'] = ids.id[ids.id.str.contains('RP')]
            
            # fill in proj_id NaNs
            ids.proj_id[ids.proj_id.isnull()] = proj_id
    
            self.log.debug(f'joining adatas for {proj_id}')
            adata = adatas[0].concatenate(adatas[1:], batch_categories = ids.id )
            adata.obs.columns = ['cell_id','id']
            ind = adata.obs.index
            tmpdf = pd.merge(adata.obs, ids, on= 'id')   # include technology in obs.
            tmpdf.index = tmpdf.cell_id
            
            adata.obs = tmpdf 
            adata.obs.index=ind
    
            # fill in missing run_ids with the cell_id. (typically smart seq experiments)
            adata.obs.run_id[adata.obs.run_id.isnull()] = adata.obs.cell_id[adata.obs.run_id.isnull()]
    
            # get batch information
            self.log.debug(f'getting batch info for {proj_id}')
            impute = load_df(f'{self.metadir}/impute.tsv')
            impute = impute.loc[impute.proj_id == proj_id,: ]
            tmpdf =  pd.merge(adata.obs, impute )[['cell_id','run_id','exp_id', 'samp_id','proj_id','tech', 'batch' ]]
            tmpdf.index= ind
            adata.obs = tmpdf
    
            adata = self._get_stats_scanpy(adata)
    
            h5file = f'{self.outputdir}/{proj_id}.h5ad'
            adata.var.index.name = None
            self.log.debug(adata.obs)
            self.log.debug(adata.var)
            self.log.debug(f'saving to {self.outputdir}/{proj_id}.h5ad ')
            adata.write_h5ad(h5file)
    
            self.log.debug(f'assigning cell types using metamarkers for {proj_id}')
            adata = self._run_MetaMarkers(h5file,adata)
            self.log.debug(f'done with metamarkers for {proj_id} - saving to h5file')
    
            adata = self._get_metadata(proj_id, adata)
            adata.write_h5ad(h5file)
            self.log.info(f'project {proj_id} done.')
            done = proj_id

        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            self.log.error(traceback.format_exc(None))
            
        return( done, part, seen)


    
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
        geneinfo = pd.read_csv(f'{self.resourcedir}/{self.species}/geneInfo.tab', sep="\t",  skiprows=1, header=None, dtype='str')
        geneinfo.columns = ['gene_accession', 'gene_symb','type_x']
        tmpdf = pd.merge(adata.var, geneinfo, on = 'gene_accession',how = 'left')        
        #tmpdf = tmpdf.fillna( 'None')
        adata.var['ribo'] = (tmpdf['type'] == "rRNA").values
        # adata.var['cytoplasm'] = None       # GO Term/kegg?
        # adata.var['metabolism'] = None      # GO Term/kegg?
        # adata.var['membrane'] = None        # GO Term/kegg?

        qcvars = ['mt', 'ribo', 'female', 'male',
                  'essential', 'cell_cycle']

        # get housekeeping genes
        self.log.debug(f"gathering marker sets:")
        hk_genes = pd.read_csv(self.stable_housekeepinig, sep=",",dtype='str')
        self.log.debug(f"housekeeping:{hk_genes}")
        adata.var['housekeeping'] = adata.var.gene_symbol.isin(
            hk_genes.gene)

        # get gender markers  https://www.nature.com/articles/s41586-020-2536-x
        # just one gene - Xist
        female_genes = pd.read_csv(self.female_genes, sep=",",dtype='str')
        self.log.debug(f"female:{female_genes}")
        adata.var['female'] = adata.var.gene_symbol.isin(
            female_genes.gene)

        # just one gene - Ddx3y
        male_genes = pd.read_csv(self.male_genes, sep=",",dtype='str')
        self.log.debug(f"male:{male_genes}")
        adata.var['male'] = adata.var.gene_symbol.isin(
            male_genes.gene)

        # get essential genes  # 1947 genes
        essential = pd.read_csv(self.essential_genes, sep=",",dtype='str')
        self.log.debug(f"ess:{essential}")
        adata.var['essential'] = adata.var.gene_symbol.isin(
            essential.gene)

        # get cell cycle genes
        cc = pd.read_csv(self.cell_cycle_genes, sep=",",dtype='str')
        cc.columns = ['cc_cluster','gene_symbol']
        self.log.debug(f"cc:{cc}")

        # get the type of gene - pseudogene, rRNA, etc
        tmpdf = pd.merge(adata.var ,cc, on = 'gene_symbol',how='left' )
        # tmpdf.drop('gene',axis=1)
        tmpdf.index = adata.var.index
        tmpdf.cc_cluster.fillna('None')
        adata.var = tmpdf
        
        
        adata.var['cell_cycle'] = adata.var.cc_cluster != 'None'
        for i in cc.cc_cluster.unique():
            adata.var[f'cc_cluster_{i}'] = adata.var.cc_cluster == i 
            qcvars.append(f'cc_cluster_{i}')

        self.log.debug('calculating qc metrics')
        sc.pp.calculate_qc_metrics(
            adata,
            expr_type='counts', var_type='genes',
            percent_top=(50, 100, 200, 500), inplace=True, use_raw=False,
            qc_vars=qcvars)

        self.log.debug('calculating corrtomean') # using log1p counts
        adata.obs['corr_to_mean'] = sparse_pairwise_corr(sparse.csr_matrix(adata.X.mean(0)) , adata.X)
        # natural log thge data


        tmp = gini_coefficient_fast(adata.X)
        adata.obs['gini']  = tmp

        sc.pp.log1p(adata)  

        # computes the N+M x N+M corrcoef matrix - extract off diagonal block
        adata.obs['log1p_corr_to_mean'] = sparse_pairwise_corr(sparse.csr_matrix(adata.X.mean(0)) , adata.X)

        # batch this.
        # X = np.random.normal(size=(10, 3))
        # F = np.zeros((10, ))

        # import multiprocessing
        # pool = multiprocessing.Pool(processes=16)
        # if number of processes is not specified, it uses the number of core
        # F[:] = pool.map(my_function, (X[i,:] for i in range(10)) )
        self.log.debug('calculating gini for each cell')
        tmp = gini_coefficient_fast(adata.X)

        adata.obs['log1p_gini']  = tmp
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
        # print('getting neighbors\n')
        sc.pp.neighbors(adata, n_neighbors=int(self.n_neighbors), n_pcs=ncomp)
        # print('getting umap\n')
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


    
    def _run_MetaMarkers(self, h5path, adata):
        '''
        h5path should be saved in the output directory before this stage. 
        After MetaMarkers, move completed h5ad file to output/
        '''


        #scriptpath = '/home/johlee/git/scqc/scqc/get_markers.R'    # should already be in path
        # should have only one outpath per dataset. 
        # TODO verify file does not already exist. Think about what to do if it does
        # outpath = f'{self.outputdir}/{os.path.basename(h5path)}'
        #cmd = ['Rscript', f'{scriptpath}',
        cmd = ['get_markers.R',
               '--class_markerset', f'{self.class_markerset}',
               '--subclass_markerset', f'{self.subclass_markerset}',
               '--h5path', f'{h5path}',
               '--max_rank', f'{self.max_rank}',
               '--outprefix', f'{self.tempdir}'
               ] 
        run_command(cmd)

        file_exts = ['_class_pred.tsv','_subclass_pred.tsv',
                     '_class_scores.tsv','_subclass_scores.tsv',
                     '_class_enrichment.tsv','_subclass_enrichment.tsv']
        tmp_paths = [ f"{self.tempdir}/{os.path.basename(h5path.replace('.h5ad',suff))}" for suff in file_exts ]
        for tmp_path in tmp_paths:
            ky = tmp_path.split('_',maxsplit=1)[1].replace('.tsv','')
            tmpdf = pd.read_csv(tmp_path ,sep="\t", index_col =0,dtype='str')
            # failsafe - some cells may be ignored in the pred file... 
            tmpdf = pd.merge(adata.obs[['cell_id']], tmpdf,how='left' , left_index = True, right_index=True)
            tmpdf = tmpdf.drop('cell_id',axis=1)
            tmpdf = tmpdf.fillna('NA')
            adata.obsm[ky] = tmpdf

            if '_class_pred.tsv' in tmp_path :
                adata.obs['class_label'] = tmpdf.predicted
            elif  '_subclass_pred.tsv' in tmp_path :
                adata.obs['subclass_label'] = tmpdf.predicted
            
            if not self.nocleanup:
                os.remove(tmp_path)
        
        return(adata)


    def _get_metadata(self, proj_id, adata):
        rdf = load_df( f'{self.metadir}/runs.tsv' )
        rdf = rdf.loc[rdf.proj_id == proj_id, :]

        edf = load_df( f'{self.metadir}/experiments.tsv' )
        edf = edf.loc[edf.proj_id == proj_id, :]

        sdf = load_df( f'{self.metadir}/samples.tsv' )
        sdf = sdf.loc[sdf.proj_id == proj_id, :]

        pdf = load_df( f'{self.metadir}/projects.tsv' )
        pdf = pdf.loc[pdf.proj_id == proj_id, :]

        adata.uns['title'] = pdf.title.values[0]
        adata.uns['abstract'] = pdf.abstract.values[0]
        adata.uns['ext_ids'] =  pdf.ext_ids.values[0]
        # TODO include external links
        # adata.uns['ext_link'] = pdf.ext_link.values[0]

        adata.uns['n_cell'] = adata.shape[0]
        adata.uns['n_runs'] = rdf.shape[0]
        adata.uns['n_exps'] = edf.shape[0]
        adata.uns['n_samp'] = sdf.shape[0]
        
        adata.uns['run_df'] = rdf
        adata.uns['exp_df'] = edf
        adata.uns['samp_df'] = sdf

        return(adata)            
    
    
    

def setup(config, overwrite=False):
    """
    create metmarker and output dirs. 
    """
    log = logging.getLogger('statistics')
    rootdir = os.path.expanduser(config.get('statistics', 'rootdir'))
    
    for sd in ['metamarker','output']:
        try:
            log.debug(f"making directory: {sd}")
            os.makedirs(f'{rootdir}/{sd}')
        except FileExistsError:
            pass
    


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
        q = Statistics(cp)
        print(args.projectids)
        for pid in args.projectids:
            # args.projectid is a list of project ids
            q.execute(pid)
