#!/usr/bin/env python
#
#  Module to deal with running STAR/StarSolo
#

import argparse
import glob
import io
import itertools
import json
import logging
import os
import re
import requests
import subprocess
import sys
import time
import urllib
import ast
from pathlib import Path
from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty
from requests.exceptions import ChunkedEncodingError

import xml.etree.ElementTree as et
import pandas as pd
import numpy as np

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import *
# inputs should be runs  identified as 'some10x'
# fastq files should already downloaded.
# srrid and species needs to be passed in from the dataframe
# can pass star parameters from config


# TODO queue 
# TODO outlists

class AlignReads(object):
    '''
    Requires:
        - Set Up to be done first.
        - fastq files for the given project in tempdir
        - STAR 2.7.# in path
        - project included in metadata
 
    '''
    def __init__(self, config ):
        self.log = logging.getLogger('star')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('star', 'tempdir'))
    
        self.metadir= os.path.expanduser(
            self.config.get('star', 'metadir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('star', 'resourcedir'))
        self.species = self.config.get('star', 'species')
        
        self.outputdir= os.path.expanduser(
            self.config.get('star', 'outputdir'))
        # self.outlist = outlist
        self.ncore_align = self.config.get('star', 'ncore_align')


    def execute(self,srpid):
        # get relevant metadata
        rdf = self._get_meta_data(srpid)
        self.log.debug(f'initializing STAR alignment for {srpid}')

        # split by technology and parses independently based on tech
        for tech, df in rdf.groupby(by = "tech_version") :
            # smartseq runs
            if tech =="smartseq":
                # build the manifest
                (manipath, manifest) = self._make_manifest(srpid,df) 
                # run star
                if manifest is not None :   # redundant ish
                    self.log.debug(f'Starting smartseq alignment for {srpid}')
                    self._run_star_smartseq(srpid,manipath)
            
            # 10x runs
            elif tech.startswith('10xv'): 

                #    run star for each run
                for row in range(df.shape[0]):
                    srrid = df.run_id[row]
                    read1 = f'{self.tempdir}/{df.read1[row]}'
                    read2 = f'{self.tempdir}/{df.read2[row]}'
                    # saves to disk
                    solooutdir = self._run_star_10x(srrid, tech, read1,read2)
                    # move solo out directories
                    
                    self._clean_up_tempdir(srpid, solooutdir )

            else :
                self.log.debug(
                    f'{tech} is not yet supported for STAR alignment.')
                # log... technology not yet supported
                pass


    # run|tech|read1|read2|exp|samp|proj|taxon|batch  dataframe in impute. 
    #       taxon to filter
    #       tech for star run
    #       reads will be used to determine which is the biological/technical (10x)
    #           read1 should be biological (cDNA)
    #           read2 should be technical (umi+cb)
    #       or for Smartseq manifest. order of reads is irrelevent.
    #       batch will be used by stats
    def _get_meta_data(self,srpid ):
        '''
        example srpid="SRP114926"
        '''

        impute = pd.read_csv(f'{self.metadir}/impute.tsv',sep="\t" ,index_col=0)

        # filter to include only requested project id
        impute = impute.loc[ impute.proj_id==srpid ,:]

        # filter to include only requested species and only keep run ids
        impute = impute.loc[ impute.taxon == int(spec_to_taxon(self.species)) ,:].reset_index(drop=True)

        return (impute)

    # smart seq scripts
    #TODO test smartseq
    def _make_manifest(self,srpid,run_data):
        # search for all fastq files with <run>_[0-9].fastq
        manipath = f"{self.tempdir}/{srpid}_smartseq_manifest.tsv"

        # runid = runlist[1]
        manifest = run_data[['read1','read2','run_id','tech_version']]
        manifest = manifest.loc[manifest.tech_version == 'smartseq',['read1','read2','run_id']]
        
        # XXX gives a warning - can't figure out how to fix....
        manifest['read1'] =  self.tempdir+'/'+manifest['read1'].astype(str)
        manifest['read2'] =  self.tempdir+'/'+manifest['read2'].astype(str)
        
        
        # overwrite
        if len(manifest) > 0 :
            manifest.to_csv(manipath, sep="\t", header=None, index=False, mode="w")
        else :
            return(manipath,None)

        return(manipath, manifest)

    def _run_star_smartseq(self,srpid,manipath):

        ss_params = {"solo_type": "SmartSeq",
                     "soloUMIdedup": "Exact",
                     "soloStrand": "Unstranded"}

        out_file_prefix = f'{self.outputdir}/{srpid}_smartseq_'
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', f'{self.ncore_align}',
               '--genomeDir', f'{self.resourcedir}/{self.species}',
               '--outFileNamePrefix', out_file_prefix,
               '--soloType', ss_params["solo_type"],
               '--soloFeatures', 'Gene',
               '--readFilesManifest', f'{manipath}',
               '--soloUMIdedup', ss_params["soloUMIdedup"],
               '--soloStrand', ss_params["soloStrand"],
               '--outSAMtype', 'None']

        cmdstr = " ".join(cmd)

        logging.debug(f"STAR command: {cmdstr} running...")
        cp = subprocess.run(cmd)

        logging.debug(
            f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")

        return( f'{out_file_prefix}Solo.out', 
                f'{out_file_prefix}Log.final.out')
    # 10x scripts
    def _get_10x_STAR_parameters(self, tech):
        d = {
            "10xv1":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv1.txt',
                    "CB_length": "14",
                    "UMI_start": "15",
                    "UMI_length": "10"
                },
            "10xv2":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv2.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "10"
                },
            "10xv3":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv3.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "12"
                }
        }
        return(d[tech])
   
    # impute stage will obtain tech, and bio/tech_readpaths for 10x runs
    def _run_star_10x(self,srrid, tech, bio_readpath,tech_readpath):
        '''
        saves to temp directory
        '''
        # read_bio, read_tech, tech = self._impute_10x_version()
        # ideally, which read is which will be obtained from impute stage
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary
        self.log.debug(f'Starting 10x alignment for {srrid}')
        cmd = ['STAR',
                '--runMode', 'alignReads',
                '--runThreadN', f'{self.ncore_align}',
                '--genomeDir', f'{self.resourcedir}/{self.species}',
                '--outFileNamePrefix', f'{self.tempdir}/{srrid}_{tech}_',
                '--soloType', star_param["solo_type"],
                '--soloCBwhitelist', star_param["white_list_path"],
                '--soloCBlen', star_param["CB_length"],
                '--soloUMIstart', f'{int(star_param["CB_length"]) + 1}',
                '--soloUMIlen', star_param["UMI_length"],
                '--soloFeatures', 'Gene',
                '--readFilesIn', bio_readpath, tech_readpath,
                '--outSAMtype', 'None']

        cp = subprocess.run(cmd)

        cmdstr = " ".join(cmd)
        logging.debug(f"STAR command: {cmdstr} running...")
        cp = subprocess.run(cmd)
        logging.debug(
            f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
        # successful runs - append to outlist.
        # if str(cp.returncode) == "0":
        #     self.outlist.append(self.srrid)

        return( f'{self.tempdir}/{srrid}_{tech}_Solo.out', 
                f'{self.tempdir}/{srrid}_{tech}_Log.final.out')

    def _clean_up_tempdir(self, srpid,solooutdir):
        try:    # make project specific solo out directories. 
            os.makedirs(f'{self.outputdir}/{srpid}')
        except FileExistsError:
            pass
        
        base = os.path.basename(solooutdir)
        dirname = os.path.dirname(solooutdir)
        newdirname = f'{self.outputdir}/{srpid}/{base}'
        try :
            os.rename(solooutdir,newdirname)
        except FileNotFoundError :
            pass

        
        starlog = base.replace('Solo.out','Log.final.out')
        newdirname = f'{self.outputdir}/{srpid}/{starlog}'
        try:
            os.rename(f'{dirname}/{starlog}',newdirname)     
        except FileNotFoundError :
            pass

        # os.remove ... Log.final.out ... Log.progress.out ...SJ.out.tab

### setup scripts
def setup(config, overwrite=False):
    '''
    Builds directories in config file 
    Download the appropriate supplement data.
    Only needs to be done once.
    '''
    log = logging.getLogger('star')
    metadir = os.path.expanduser(config.get('star', 'metadir'))
    cachedir = os.path.expanduser(config.get('star', 'cachedir'))
    tempdir = os.path.expanduser(config.get('star', 'tempdir'))
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    #n_core = config.get('star', 'n_core')

    for d in [metadir, cachedir, tempdir, resourcedir]:
        try:
            log.debug(f"making directory: {d}")
            os.makedirs(d)
        except FileExistsError:
            pass

    get_whitelists(config, overwrite)
    get_genome_data(config, overwrite)
    build_genome_indices(config, overwrite)


def get_whitelists(config,force=False):
    """
    Get cellranger tag whitelists. Assumes resourcedir exists. 
    """
    log = logging.getLogger('star')
    wls = ['whitelist_10xv1', 'whitelist_10xv2', 'whitelist_10xv3']
    outdir = os.path.expanduser(config.get('star', 'resourcedir'))

    for key in wls:
        url = config.get('star', key)
        log.debug(f'getting cellranger whitelist from {url}...')
        r = requests.get(url)
        if r.status_code == 200:
            with open("/".join([outdir, key]), "w") as f:
                if url.endswith(".gz"):
                    f.write(gzip.decompress(r.content).decode())
                else:
                    f.write(r.text)
        else:
            log.warning(
                f"Retrieving {url} failed with status_code: {str(r.status_code)}")


def get_genome_data(config, overwrite=False):
    """
    Download required genome data
    Assumes config is for one species only. 
    
    """
    log = logging.getLogger('star')
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    species = config.get('star', 'species')
    log.debug(f'got species {species}')
    outdir = "/".join([resourcedir, species])
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    fa_url = config.get('star', f'{species}_fa_url')
    log.debug(f'{species} fa_url={fa_url}')
    download_ftpurl(fa_url, outdir, 'genome.fa', overwrite)

    gtf_url = config.get('star', f'{species}_gtf_url')
    log.debug(f'{species} gtf_url={gtf_url}')
    download_ftpurl(gtf_url, outdir, 'annotation.gtf', overwrite)


def build_genome_indices(config, force=False):
    log = logging.getLogger('star')

    n_core = int(config.get('star', 'ncore_index'))
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    speciesnames = config.get('star', 'species')
    specieslist = [i.strip() for i in speciesnames.split(',')]
    log.debug(f'got species {specieslist}')

    for species in specieslist:
        outdir = "/".join([resourcedir, species])
        # NOTE do we care about low memory here? Only done once and we have the RAM
        cmd = ["STAR",
               "--runMode", "genomeGenerate",
               "--genomeSAsparseD", "3",   # for low memory (RAM)
               "--genomeSAindexNbases", "12"  # for low memory
               "--runThreadN", f'{n_core}',
               "--genomeDir", f'{outdir}',
               "--genomeFastaFiles", f'{outdir}/genome.fa',
               "--sjdbGTFfile", f'{outdir}/annotation.gtf']

        cmdstr = " ".join(cmd)
        log.info(f'running {cmdstr} ...')
        cp = subprocess.run(cmdstr,
                            shell=True,
                            universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        log.debug(f'Ran cmd={cmdstr}  returncode={cp.returncode} ')

    logging.debug(
        f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
    
    if str(cp.returncode) != "0":
        log.warning(f"non-zero return code from STAR. See Log.out...")



# to do: include drivers for 10x and ss alignments
if __name__ == "__main__":

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
 
    parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default='~/git/scqc/etc/scqc.conf',
                            help='Config file path [~/git/scqc/etc/scqc.conf]')

    parser.add_argument('-s', '--setup',
                        action='store_true',
                        dest='setup',
                        help='Set up directories and downloads supplemental data')

    parser.add_argument('-F', '--force',
                        action='store_true',
                        dest='force',
                        default=False,
                        help='Re-do set up, overwriting what is done.')

    parser.add_argument('-p', '--project_id',
                        metavar='project_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Align reads for all runs in project_id')

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
    
    logging.debug(f"args: {args}")
    logging.debug(f"got config: {cs}")

    if args.setup:
        s = setup(cp, overwrite=args.force)
        # s.execute()

    
    if args.project_id:
        a = AlignReads(cp)
        for pid in args.project_id:
            a.execute(pid)
