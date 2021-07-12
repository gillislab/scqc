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
# import matplotlib.cbook as cbook

plt.rcParams.update({'figure.autolayout': True})

sns.set_style("ticks")
sns.set_context("paper")
sns.despine()

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

class BuildFigures(object):
    def __init__(self,config):
        self.log = logging.getLogger('figures')
        self.config = config
        self.figuredir = os.path.expanduser(    # where to store figures
            self.config.get('figures', 'figuredir'))
        self.outputdir = os.path.expanduser(    # where the adata is stored
            self.config.get('figures', 'outputdir'))

    def execute(self, srpid):
        # get the srpid

        # get the h5ad paths for the project
        
        h5paths = glob.glob(f'{self.outputdir}/{srpid}.h5ad')

        try:
            os.makedirs(f'{self.figuredir}/{srpid}')
        except FileExistsError:
            pass

        for h5path in h5paths:
            basename = os.path.basename(h5path)
            rootname = os.path.dirname(h5path)

            adata = sc.read_h5ad(h5path)

            # plot var x against var y
            figurepath = basename.replace('.h5ad', '_Gender.png')
            plot_scatter(data=adata.obs, x_column='pct_counts_male', y_column='pct_counts_female',  figsize=(4,3),
                        xlab='Pct Male Markers', ylab=' Pct Male Markers', 
                        color_by='pct_counts_male',cbar_title='Pct Male Markers',
                        regression_line=False, hlines = [], vlines =[], 
                        save_it=True, outfile=f'{self.figuredir}/{srpid}/{figurepath}')
            
            # plot on umap 
            umapdf = pd.DataFrame({
                    'x': adata.obsm['X_umap'][:,0],
                    'y':adata.obsm['X_umap'][:,1], 
                    'Corr_to_Mean':adata.obs['corr_to_mean'],
                    'Mito':adata.obs['log1p_total_counts_mt'],
                    # 'Batch':adata.obs['batch'],
                    'Cell_Cycle':adata.obs['log1p_total_counts_cell_cycle'],
                    'Essential':adata.obs['log1p_total_counts_essential'],
                    'Ribosome':adata.obs['log1p_total_counts_ribo'],
                    'Run_ID':adata.obs['run_id'],
                    'Experiment_ID':adata.obs['exp_id'],
                    'Sample_ID':adata.obs['samp_id'],
                    'Cell_ID':adata.obs['cell_id']
                    })

            for kw in umapdf.columns[2:] :
                figurepath = basename.replace('.h5ad', f'_umap_{kw}.png')
                plot_scatter(data=umapdf, x_column='x', y_column='y',  figsize=(4,3),
                            xlab='UMAP_1', ylab=' UMAP_2', 
                            color_by=kw,cbar_title=kw,
                            regression_line=False, hlines = [], vlines =[], 
                            save_it=True, outfile=f'{self.figuredir}/{srpid}/{figurepath}')

            # plot pairplots 
            cols= [ 'log1p_n_genes_by_counts','log1p_total_counts', 'log1p_total_counts_mt',
                    'log1p_total_counts_essential', 
                    'log1p_total_counts_cell_cycle','corr_to_mean', 'gini']
            sns.pairplot(adata.obs[cols])
            plt.savefig(f'{self.figuredir}/{srpid}_pairplot.png'  , dpi=300, transparent=False )


            # plot histograms 
            figurepath = basename.replace('.h5ad',f'_gini.png')
            plot_histogram(adata.obs.gini    , xlab='Gini coefficient', ylab='Frequency', title=None, save_it=True, outfile= f'{self.figuredir}/{srpid}/{figurepath}',
                  bins=None, orientation='vertical', vlines=['mean'], figsize=(4, 3), rotate_xlabel=False)
            figurepath = basename.replace('.h5ad',f'_Corr_to_Mean.png')
            plot_histogram(adata.obs.corr_to_mean    , xlab='Correlation to Mean', ylab='Frequency', title=None, save_it=True, outfile= f'{self.figuredir}/{srpid}/{figurepath}',
                  bins=None, orientation='vertical', vlines=['mean'], figsize=(4, 3), rotate_xlabel=False)


# sns.pairplot(adata.obs[['log1p_n_genes_by_counts', 'log1p_total_counts',
#              'log1p_total_counts_mt',
#              'log1p_total_counts_ribo', 'log1p_total_counts_essential',
#              'log1p_total_counts_cell_cycle']] , height = 6, aspect=1)

def plot_scatter(data, x_column, y_column, xlab=None, ylab=None, title =None, figsize=(4,3), 
                marker='.',markersize = 1, color_by=None, cbar_title=None,
                regression_line=True,  hlines = [], vlines =[],
                save_it=False, outfile="/home/johlee/scqc/figures/temp_sns.png"):

    colors = None
    cmap =None
    fig, ax = plt.subplots(figsize=figsize)


    if color_by  is not None :
        print(color_by)
        colors = list(data[color_by])
        cmap = sns.color_palette("mako_r", as_cmap=True)    # for continuous data


        ty = type(colors[0])
        if  ty == str:
            if len(set(colors)) ==  len(colors) or len(set(colors))> 20 :    # categorical - all unique or too many cats
                scat = ax.scatter(x=data[x_column], y=data[y_column],  marker=marker,s = markersize )
            else :  
                for name, group in data.groupby(color_by):
                    print(name)
                    scat = ax.scatter(x=group[x_column], y=group[y_column],  marker=marker,s = markersize, label=name )
                    ax.legend(bbox_to_anchor=(1, 1),fontsize='xx-small',title=color_by)

        else :
            scat = ax.scatter(x=data[x_column], y=data[y_column],  c= colors, marker=marker,s = markersize, cmap = cmap)

            cb = plt.colorbar(scat)
            cb.set_label(cbar_title)
    
    scat = ax.set(xlabel=xlab, ylabel=ylab,
           title=title)

    if regression_line :    # add a regression line
        m, b = np.polyfit(data[x_column], data[y_column], 1)
        plt.plot(data[x_column],data[x_column]*m+b, 'k--' ,linewidth=markersize-.5)
        # sns.regplot(,y,ci=None)

    for vline in list(vlines):
        if vline == 'mean':
            m = data[x_column].mean()
            ax.text(m, 0, 'Mean', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r',linewidth=markersize-.5)
        elif vline == 'median':
            m = data[x_column].median()
            ax.text(m, 0, 'Median', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r',linewidth=markersize-.5)
        elif type(vline) == int or type(vline) == float:
            ax.axvline(float(vline), ls='--', color='r',linewidth=markersize-.5)


    for hline in list(hlines):
        if hline == 'mean':
            m = data[y_column].mean()
            ax.text(m, 0, 'Mean', fontsize=10)
            ax.axyline(m, ls='--', color='r',linewidth=markersize-.5)
        elif vline == 'median':
            m = data[y_column].median()
            ax.text(m, 0, 'Median', fontsize=10)
            ax.axhline(m, ls='--', color='r',linewidth=markersize-.5)
        elif type(hline) == int or type(vline) == float:
            ax.axhline(float(vline), ls='--', color='r',linewidth=markersize-.5)

    if save_it:
        # transparent should be set to true for publication etc. False for testing.
        fig.savefig(outfile, dpi=300, transparent=False)

    return( fig, ax )

def plot_histogram(x, xlab=None, ylab=None, title=None, save_it=False, outfile='/home/johlee/scqc/figures/temp_hist.png',
                  bins=None, orientation='vertical', vlines=['mean', 'median'], figsize=(8, 6), rotate_xlabel=False ):
    fig, ax = plt.subplots(figsize=figsize)

    ax.hist(x, bins=bins, orientation=orientation)
    ax.set(xlabel=xlab, ylabel=ylab,
           title=title)

    for vline in list(vlines):

        if vline == 'mean':
            m = x.mean()
            ax.text(m, 0, 'Mean', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r')
        elif vline == 'median':
            m = x.median()
            ax.text(m, 0, 'Median', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r')
        elif type(vline) == int or type(vline) == float:
            ax.axvline(float(vline), ls='--', color='r')

    for vline in list(vlines):

        if vline == 'mean':
            m = x.mean()
            ax.text(m, 0, 'Mean', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r')
        elif vline == 'median':
            m = x.median()
            ax.text(m, 0, 'Median', fontsize=10, horizontalalignment='center')
            ax.axvline(m, ls='--', color='r')
        elif type(vline) == int or type(vline) == float:
            ax.axvline(float(vline), ls='--', color='r')

    if rotate_xlabel:
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=45, horizontalalignment='center')

    if save_it:
        # transparent should be set to true for publication etc. False for testing.
        fig.savefig(outfile, dpi=300, transparent=False)

    return(fig, ax)

def plot_boxplots(data): 
    # primarily for cell cycle phases
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
        q = BuildFigures(cp)
        for pid in args.projectids:
            # args.projectid is a list of project ids
            q.execute(pid)
