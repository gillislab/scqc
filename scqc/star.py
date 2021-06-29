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



class AlignReads(object):
    '''
    Requires:
        - Set Up to be done first.
        - fastq files for the given project in tempdir
        - STAR 2.7.# in path
        - project included in metadata
 
    '''
    def __init__(self, config, srpid ):
        self.log = logging.getLogger('star')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('star', 'tempdir'))
    
        self.metadir= os.path.expanduser(
            self.config.get('star', 'metadir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('star', 'resourcedir'))
        self.species = self.config.get('star', 'tempdir')
        
        self.outputdir= os.path.expanduser(
            self.config.get('star', 'outputdir'))
        # self.outlist = outlist
        self.ncore_align = self.config.get('star', 'ncore_align')

        self.log.debug(f'initializing STAR alignment for {srpid}')
        self.srpid=srpid


    def execute(self):
        # get relevant metadata
        rdf = self._get_meta_data()
        # split by technology and parses independently based on tech
        for tech, df in rdf.groupby(by = "tech") :
            # smartseq runs
            if tech =="smartseq":
                # build the manifest
                (manipath, manifest) = self._make_manifest(df) 
                # run star
                self._run_star_smartseq(manipath,manifest)
            
            # 10x runs
            if tech.startswith('10xv'): 
                # grab star params
                # for run in run_ids
                #    run star
                for row in df.shape[0]:
                    srrid = df.run_id[row]
                    read1 = df.read1[row]
                    read2 = df.read2[row]
                    # saves to disk
                    self._run_star_10x(srrid, tech, read1,read2)
                pass 
            else :
                self.log.debug(
                    f'{tech} is not yet supported for STAR alignment.')
                # log... technology not yet supported
                pass


    # TODO  make a run|taxon|tech|read1|read2|batch  dataframe in impute. 
    #       taxon to filter
    #       tech for star run
    #       reads will be used to determine which is the biological/technical
    #           read1 should be biological (cDNA)
    #           read2 should be technical (umi+cb)
    #       batch will be used by stats
    def _get_meta_data(self ):
        '''
        example srpid="SRP114926"
        '''
        
        # sdf = pd.read_csv(f'{metadir}/samples.tsv',sep="\t" ,index_col=0)
        # edf = pd.read_csv(f'{metadir}/experiments.tsv',sep="\t" ,index_col=0)
        rdf = pd.read_csv(f'{self.metadir}/runs.tsv',sep="\t" ,index_col=0)
        
        # TODO from imputation - df containing runs with corresponding tech
        run2tech = pd.read_csv(f'{self.metadir}/run2tech.tsv',sep="\t" ,index_col=0)
        run2tech.columns = ["run_id","tech"]
        # filter to include only requested project id
        rdf = rdf.loc[ rdf.proj_id==self.srpid ,:]
        # filter to include only requested species and only keep run ids
        rdf = rdf.loc[ rdf.taxon == int(spec_to_taxon(self.species)) ,['run_id','nreads']]

        rdf = pd.merge(rdf, run2tech , how = 'left', on = "run_id") 

        return (rdf)

    # smart seq scripts
    def _make_manifest(self,run_data):
        # search for all fastq files with <run>_[0-9].fastq
        manipath = f"{self.metadir}/{self.srpid}_smartseq_manifest.tsv"

        # runid = runlist[1]
        allRows = []
        for runid in run_data.run_id:
            # where are the fastq files? In temp
            fqs = glob.glob(f'{self.tempdir}/{runid}*.fastq')
            fqs.sort()
            # number of fastq files found for the run
            if len(fqs) > 0 and len(fqs) < 3:
                if len(fqs) == 1:
                    fqs.append('-')
                    fqs.append(runid)
                elif len(fqs) == 2:
                    fqs.append(runid)

                allRows.append(fqs)

        manifest = pd.DataFrame(allRows, columns=['read1', 'read2', 'run'])

        # overwrite
        manifest.to_csv(manipath, sep="\t", header=None, index=False, mode="w")

        return(manipath, manifest)

    def _run_star_smartseq(self,manifest,manipath):

        ss_params = {"solo_type": "SmartSeq",
                     "soloUMIdedup": "Exact",
                     "soloStrand": "Unstranded"}

        # TODO filter to the runs that we havne't aligned yet. 
        #   adjust manifest accordingly.

        # runs_done_file = f'{self.outputdir}/{self.srpid}_smartseq_Solo.out/Gene/raw/barcodes.tsv'
        # if os.path.isfile(runs_done_file):
        #     runs_done = open(runs_done_file).read().strip().split('\n')

        #     # of the runs that i find in the manifest, which have already been aligned?
        #     new_runs = listdiff(manifest.run.values, runs_done)
        #     #list(set(manifest.run.values) - set(runs_done)).sort()
        #     tmp_mani = manifest.loc[new_runs == manifest.run, :]

        #     # save tmp manifest to temp directory - may be empty
        #     tmp_manipath = manipath.replace(
        #         f'{self.metadir}', f'{self.tempdir}')
        #     tmp_mani.to_csv(tmp_manipath, sep="\t")
        #     out_file_prefix = f'{self.tempdir}/{self.srpid}_smartseq_'
        # else:
        #     tmp_manipath = manipath
        #     out_file_prefix = f'{self.staroutdir}/{self.srpid}_smartseq_'
        
        out_file_prefix = f'{self.outputdir}/{self.srpid}_smartseq_'
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', f'{self.num_streams}',
               '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
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
        # successful runs - append to outlist.
        if str(cp.returncode) == "0":
            self.outlist.append(self.srrid)

        # # did we write to a temp directory?
        # if out_file_prefix.startswith(f'{self.tempdir}'):
        #     self._merge_solo_out_results(
        #         f'{self.staroutdir}/{self.srpid}_smartseq_',    # starout direc
        #         out_file_prefix)                                # temp direc

    # 10x scripts
    def _get_10x_STAR_parameters(self, tech):
        d = {
            "10xv1":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv1.txt',
                    "CB_length": "14",
                    "UMI_start": "15",
                    "UMI_length": "10"
                },
            "10xv2":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv2.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "10"
                },
            "10xv3":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv3.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "12"
                }
        }
        return(d[tech])
   
    # impute stage will obtain tech, and bio/tech_readpaths for 10x runs
    def _run_star_10x(self,srrid, tech, bio_readpath,tech_readpath):

        # read_bio, read_tech, tech = self._impute_10x_version()
        # ideally, which read is which will be obtained from impute stage
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary

        cmd = ['STAR',
                '--runMode', 'alignReads',
                '--runThreadN', f'{self.ncore_align}',
                '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
                '--outFileNamePrefix', f'{self.outputdir}/{srrid}_{tech}_',
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
        logging.debug(f"Fasterq-dump command: {cmdstr} running...")
        cp = subprocess.run(cmd)
        logging.debug(
            f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
        # successful runs - append to outlist.
        if str(cp.returncode) == "0":
            self.outlist.append(self.srrid)

 




class Align10xSTAR(object):
    '''
        - Identifies 10x version
        - gets the path to fastq files
        - 
        Simple wrapper for STAR - 10x input
    '''

    def __init__(self, config, srrid, species, outlist):
        self.log = logging.getLogger('star')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('analysis', 'tempdir'))
        self.srrid = srrid
        self.log.debug(f'aligning id {srrid}')
        self.staroutdir = os.path.expanduser(
            self.config.get('analysis', 'staroutdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('analysis', 'resourcedir'))
        self.species = species
        self.outlist = outlist
        self.num_streams = self.config.get('analysis', 'num_streams')


    # tested on SRR14633482 - did not get a Solo.out directory? Ran with 10xv3 params (though umi+cb=30)
    def execute(self):
        read_bio, read_tech, tech = self._impute_10x_version()
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary

        if tech != 'unknown':
            cmd = ['STAR',
                   '--runMode', 'alignReads',
                   '--runThreadN', f'{self.num_streams}',
                   '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
                   '--outFileNamePrefix', f'{self.staroutdir}/{self.srrid}_{tech}_',
                   '--soloType', star_param["solo_type"],
                   '--soloCBwhitelist', star_param["white_list_path"],
                   '--soloCBlen', star_param["CB_length"],
                   '--soloUMIstart', f'{int(star_param["CB_length"]) + 1}',
                   '--soloUMIlen', star_param["UMI_length"],
                   '--soloFeatures', 'Gene',
                   '--readFilesIn', read_bio, read_tech,
                   '--outSAMtype', 'None']

            cp = subprocess.run(cmd)

            cmdstr = " ".join(cmd)
            logging.debug(f"Fasterq-dump command: {cmdstr} running...")
            cp = subprocess.run(cmd)
            logging.debug(
                f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
            # successful runs - append to outlist.
            if str(cp.returncode) == "0":
                self.outlist.append(self.srrid)

        else:
            pass


# currently, if solo.out dir exsts, save results to the temp directory
# then merge results with previous star runs (to do)
# input should be a project accession
# no real way to verify that the reads are indeed smartseq...
# can pass star parameters from config

def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp


# john lee is satisfied with this class 6/3/2021
def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()



class AlignSmartSeqSTAR(object):

    def __init__(self, config, species, srpid, outlist):
        self.log = logging.getLogger('sra')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('analysis', 'tempdir'))
        self.srpid = srpid
        self.log.debug(f'aligning smartseq run from {srpid}')
        self.staroutdir = os.path.expanduser(
            self.config.get('analysis', 'staroutdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('analysis', 'resourcedir'))
        self.species = species
        self.outlist = outlist
        self.num_streams = self.config.get('analysis', 'num_streams')

    def _make_manifest(self):
        # search for all fastq files with <run>_[0-9].fastq
        manipath = f"{self.metadir}/{self.srpid}_smartseq_manifest.tsv"

        # metadata here should only contain smart seq data.
        df = pd.read_csv(
            f'{self.metadir}/{self.srpid}_smartseq_metadata.tsv', sep="\t")

        df.runs = df.runs.apply(ast.literal_eval)
        df = df.explode('runs')
        # runlist = df.runs.values

        # runid = runlist[1]
        allRows = []
        for runid in df.runs:
            fqs = glob.glob(f'{self.cachedir}/{runid}*.fastq').sort()

            if len(fqs) > 0 and len(fqs) < 3:  # fastq files found
                if len(fqs) == 1:
                    fqs.append(['-', runid])
                elif len(fqs) == 2:
                    fqs.append(runid)

                allRows.append(fqs)

        manifest = pd.DataFrame(allRows, columns=['read1', 'read2', 'run'])

        # overwrite
        manifest.to_csv(manipath, sep="\t", header=None, index=False, mode="w")

        return(manipath, manifest)

    # to do;
    def _merge_solo_out_results(self, final_path, temp_path):
        '''
        Merge the Solo.out results for two star runs on different parts of the data
        '''
        # merge data,

        # delete from temp

        pass

    def execute(self):

        ss_params = {"solo_type": "SmartSeq",
                     "soloUMIdedup": "Exact",
                     "soloStrand": "Unstranded"}

        manipath, manifest = self._make_manifest()

        # do we already have star output for this project.
        # If so, save starout results in a temp directory, merge data/stats
        # then delete tmp
        runs_done_file = f'{self.staroutdir}/{self.srpid}_smartseq_Solo.out/Gene/raw/barcodes.tsv'
        if os.path.isfile(runs_done_file):
            runs_done = open(runs_done_file).read().strip().split('\n')

            # of the runs that i find in the manifest, which have already been aligned?
            new_runs = listdiff(manifest.run.values, runs_done)
            #list(set(manifest.run.values) - set(runs_done)).sort()
            tmp_mani = manifest.loc[new_runs == manifest.run, :]

            # save tmp manifest to temp directory - may be empty
            tmp_manipath = manipath.replace(
                f'{self.metadir}', f'{self.tempdir}')
            tmp_mani.to_csv(tmp_manipath, sep="\t")
            out_file_prefix = f'{self.tempdir}/{self.srpid}_smartseq_'
        else:
            tmp_manipath = manipath
            out_file_prefix = f'{self.staroutdir}/{self.srpid}_smartseq_'

        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', f'{self.num_streams}',
               '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
               '--outFileNamePrefix', out_file_prefix,
               '--soloType', ss_params["solo_type"],
               '--soloFeatures', 'Gene',
               '--readFilesManifest', f'{tmp_manipath}',
               '--soloUMIdedup', ss_params["soloUMIdedup"],
               '--soloStrand', ss_params["soloStrand"],
               '--outSAMtype', 'None']

        cmdstr = " ".join(cmd)

        logging.debug(f"STAR command: {cmdstr} running...")
        cp = subprocess.run(cmd)

        logging.debug(
            f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
        # successful runs - append to outlist.
        if str(cp.returncode) == "0":
            self.outlist.append(self.srrid)

        # did we write to a temp directory?
        if out_file_prefix.startswith(f'{self.tempdir}'):
            self._merge_solo_out_results(
                f'{self.staroutdir}/{self.srpid}_smartseq_',    # starout direc
                out_file_prefix)                                # temp direc


### setup scripts
def setup(config, force=False):
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

    get_whitelists(config, force)
    get_genome_data(config, force)
    build_genome_indices(config, force)


def get_whitelists(config,force=False):
    """
    Get cellranger tag whitelists. Assumes resourcedir exists. 
    """
    log = logging.getLogger('star')
    wls = ['10x_v1_whitelist', '10x_v2_whitelist', '10x_v3_whitelist']
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


def get_genome_data(config,force=False):
    """
    Download required genome data
    """
    log = logging.getLogger('star')
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    speciesnames = config.get('star', 'species')
    specieslist = [i.strip() for i in speciesnames.split(',')]
    log.debug(f'got species {specieslist}')

    for species in specieslist:
        outdir = "/".join([resourcedir, species])
        try:
            os.makedirs(outdir)
        except FileExistsError:
            pass

        fa_url = config.get('star', f'{species}_fa')
        log.debug(f'{species} fa_url={fa_url}')
        download_ftpurl(fa_url, outdir, 'genome.fa')

        gtf_url = config.get('star', f'{species}_gtf')
        log.debug(f'{species} gtf_url={gtf_url}')
        download_ftpurl(gtf_url, outdir, 'annotation.gtf')


def build_genome_indices(config, force=False):
    log = logging.getLogger('star')

    n_core = int(config.get('star', 'ncore_index'))
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    speciesnames = config.get('star', 'species')
    specieslist = [i.strip() for i in speciesnames.split(',')]
    log.debug(f'got species {specieslist}')

    for species in specieslist:
        outdir = "/".join([resourcedir, species])
        cmd = ["STAR",
               "--runMode", "genomeGenerate",
               "--genomeSAsparseD", "3",   # for low memory
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
    if str(cp.returncode) == "0":
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

    parser.add_argument('-s', '--setup',
                        action='store_true',
                        dest='setup',
                        help='Set up directories and downloads supplemental data')

    parser.add_argument('-F', '--force',
                        action='store_true',
                        dest='force',
                        default=False,
                        help='Re-do, overwriting what is done.')

    parser.add_argument('-f', '--fasterq',
                        metavar='fasterq',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')

    parser.add_argument('-tx', '--tenx',
                        metavar='tenx_align',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Align 10x args with STAR. e.g. SRR14584407')

    parser.add_argument('-ss', '--smartseq',
                        metavar='ss_align',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Align SmartSeq args with STAR. e.g. SRP308826')

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    cp = get_default_config()
    cs = get_configstr(cp)

    logging.debug(f"got config: {cs}")

    if args.setup:
        s = setup(cp, force=args.force)
        s.execute()
