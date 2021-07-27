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
import shutil
import subprocess
import sys
import time
import urllib
import ast
import datetime as dt
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
from scqc.sra import FasterqDump

# inputs should be runs  identified as 'some10x'
# .sra files should already downloaded.
# run_id and species needs to be passed in from the dataframe
# can pass star parameters from config


# TODO queue 

class UnsupportedTechnologyException(Exception):
    """
    Thrown when run technology is neither 10x nor SmartSeq
    """
class FasterqFailureException(Exception):
    """
    Thrown when run technology is neither 10x nor SmartSeq
    """



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
        self.cachedir = os.path.expanduser(
            self.config.get('star', 'cachedir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('star', 'resourcedir'))
        self.species = self.config.get('star', 'species')
        self.ncore_align = self.config.get('star', 'ncore_align')



    def execute(self, proj_id):
        # get relevant metadata
        rdf = self._get_meta_data(proj_id)
        self.log.debug(f'initializing STAR alignment for {proj_id}')
        done = None
        seen = proj_id
        try:
            # bring in all fastqs to <tempdir>
            runlist = self._stage_in(proj_id, rdf)
            # split by technology and parses independently based on tech
            for tech, df in rdf.groupby(by = "tech_version") :
                # smartseq runs
                if tech =="smartseq":
                    self._handle_smartseq(proj_id, df)
                elif tech.startswith('10xv'):     
                    self._handle_10x(proj_id, tech, df)
                else :
                    self.log.warning(
                        f'{tech} is not yet supported for STAR alignment.')
                    # log... technology not yet supported
                    raise UnsupportedTechnologyException(f'For project {proj_id}')                
                # finally, clean up fastq files             
            self._remove_fastqs(runlist)
            done = proj_id
        
        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            self.log.error(traceback.format_exc(None))
            
        finally:
            return (done, seen)



    def _stage_in(self, proj_id, rdf):
        """
        bring in fastq files to <tempdir> for all exp_ids in this project.
        
        throws FasterqFailureException if there is a problem. 
        
        """
        runlist  = list( rdf[ rdf.proj_id==proj_id ].run_id.unique()  ) 
        runlength = len(runlist)
        i = 0
        for run_id in runlist:
            i += 1
            fqd = FasterqDump(self.config, run_id)
            rc = fqd.execute()
            if str(rc)!= '0':
                raise FasterqFailureException(f'runid {run_id}')
            else:
                self.log.info(f'runid {run_id}  [{i}/{runlength}] handled successfully.')
        self.log.info(f'successfully extracted all {runlength} fastqs for project {proj_id}.') 
        return runlist
        


    def _handle_smartseq(self, proj_id, df):
        # build the manifest
        (manipath, manifest) = self._make_manifest(proj_id, df) 
        # run star
        self.log.debug(f'Starting smartseq alignment for {proj_id}')
        outfile_prefix = self._run_star_smartseq(proj_id, manipath)
        self.log.debug(f'Got outfile_prefix={outfile_prefix} for {proj_id}')
        self._stage_out(proj_id, outfile_prefix)
        


    def _handle_10x(self, proj_id, tech, df):
        for row in range(df.shape[0]):
            run_id = df.run_id[row]
            read1 = f'{self.tempdir}/{df.read1[row]}' # biological cDNA
            read2 = f'{self.tempdir}/{df.read2[row]}' # technical CBarcode + UMI
            # saves to disk
            outfile_prefix = self._run_star_10x(run_id, tech, read1, read2)
            self.log.debug(f'Got outfile_prefix={outfile_prefix} for {proj_id} and {run_id}')
            self._stage_out(proj_id, outfile_prefix)
            
        

    # run|tech|read1|read2|exp|samp|proj|taxon|batch  dataframe in impute. 
    #       taxon to filter
    #       tech for star run
    #       reads will be used to determine which is the biological/technical (10x)
    #           read1 should be biological (cDNA)
    #           read2 should be technical (umi+cb)
    #       or for Smartseq manifest. order of reads is irrelevent.
    #       batch will be used by stats
    def _get_meta_data(self, proj_id ):
        '''
        example proj_id="SRP114926"
        '''
        impute = pd.read_csv(f'{self.metadir}/impute.tsv',sep="\t" ,index_col=0)
        # filter to include only requested project id
        impute = impute.loc[ impute.proj_id==proj_id ,:]
        # filter to include only requested species and only keep run ids
        impute = impute.loc[ impute.taxon == int(spec_to_taxon(self.species)) ,:].reset_index(drop=True)
        return(impute)


    # smart seq scripts
    def _make_manifest(self, proj_id, run_data):
        # search for all fastq files with <run>_[0-9].fastq
        manipath = f"{self.tempdir}/{proj_id}_smartseq_manifest.tsv"

        # runid = runlist[1]
        manifest = run_data[['read1','read2','run_id','tech_version']]
        manifest = manifest.loc[manifest.tech_version == 'smartseq',['read1','read2','run_id']]       
        manifest['read1'] =  self.tempdir+'/'+manifest['read1'].astype(str)
        manifest['read2'] =  self.tempdir+'/'+manifest['read2'].astype(str)
        
        # overwrite
        if len(manifest) > 0 :
            manifest.to_csv(manipath, sep="\t", header=None, index=False, mode="w")
        else :
            return(manipath,None)

        return(manipath, manifest)


    def _run_star_smartseq(self, proj_id, manipath):

        ss_params = {"solo_type": "SmartSeq",
                     "soloUMIdedup": "Exact",
                     "soloStrand": "Unstranded"}

        outfile_prefix = f'{self.tempdir}/{proj_id}_smartseq_'
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', f'{self.ncore_align}',
               '--genomeDir', f'{self.resourcedir}/{self.species}',
               '--outFileNamePrefix', outfile_prefix,
               '--soloType', ss_params["solo_type"],
               '--soloFeatures', 'Gene',
               '--readFilesManifest', f'{manipath}',
               '--soloUMIdedup', ss_params["soloUMIdedup"],
               '--soloStrand', ss_params["soloStrand"],
               '--outSAMtype', 'None']

        self._run_command(cmd)
        return(outfile_prefix)

    
    # 10x scripts
    def _get_10x_STAR_parameters(self, tech):
        d = {
            "10xv1":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv1',
                    "CB_length": "14",
                    "UMI_start": "15",
                    "UMI_length": "10"
                },
            "10xv2":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv2',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "10"
                },
            "10xv3":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelist_10xv3',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "12"
                }
        }
        return(d[tech])
   
    # impute stage will obtain tech, and bio/tech_readpaths for 10x runs
    def _run_star_10x(self, run_id, tech, bio_readpath, tech_readpath):
        '''
        saves to temp directory
        '''
        # read_bio, read_tech, tech = self._impute_10x_version()
        # ideally, which read is which will be obtained from impute stage
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary
        self.log.debug(f'Starting 10x alignment for {run_id}')
        outfile_prefix = f'{self.tempdir}/{run_id}_{tech}_'
        cmd = ['STAR',
                '--runMode', 'alignReads',
                '--runThreadN', f'{self.ncore_align}',
                '--genomeDir', f'{self.resourcedir}/{self.species}',
                '--outFileNamePrefix', outfile_prefix ,
                '--soloType', star_param["solo_type"],
                '--soloCBwhitelist', star_param["white_list_path"],
                '--soloCBlen', star_param["CB_length"],
                '--soloUMIstart', f'{int(star_param["CB_length"]) + 1}',
                '--soloUMIlen', star_param["UMI_length"],
                '--soloFeatures', 'Gene',
                '--readFilesIn', bio_readpath, tech_readpath,
                '--outSAMtype', 'None']
        self._run_command(cmd)
        return(outfile_prefix)


    def _run_command(self, cmd):
        cmdstr = " ".join(cmd)
        self.log.info(f"command: {cmdstr} running...")
        start = dt.datetime.now()
        cp = subprocess.run(cmd, 
                        text=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT)
        end = dt.datetime.now()
        elapsed =  end - start
        self.log.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
        
        if str(cp.returncode) == '0':
            self.log.info(f'successfully ran {cmdstr}')
        else:
            self.log.error(f'non-zero return code for {self.run_id}')
        if cp.stderr is not None:
            self.log.debug(f"got stderr: {cp.stderr}")
        if cp.stdout is not None:
            self.log.debug(f"got stdout: {cp.stdout}")
            
        

    def _stage_out(self, proj_id, outfile_prefix):
        """
        Stage out. 
        <tempdir>/{outfile_prefix}Solo.out
                  {outfile_prefix}Log.out
                  {outfile_prefix}Log.final.out
                  {outfile_prefix}SJ.out.tab
                  {outfile_prefix}manifest.tsv        
        
        for 10x
        f'{self.tempdir}/{run_id}_{tech}_Solo.out' -> SRR10285015_10xv2_Solo.out
            -> <cachedir>/proj_id/{run_id}_{tech}_Solo.out/
            
        for smartseq
        f'{self.tempdir}/{proj_id}_smartseq_Solo.out' -> SRP066963_smartseq_Solo.out/
            -> <cachedir>/proj_id/{proj_id}_{tech}_Solo.out/
        
        """
        
        MOVES = ['Log.out','Log.final.out','SJ.out.tab','manifest.tsv']
          
        projdir = f'{self.cachedir}/{proj_id}/'
        
        try:    # make project specific solo out directories. 
            os.makedirs(projdir)
            self.log.debug(f'created cache project dir {projdir}')
        except FileExistsError:
            self.log.warning(f'cache project dir already exists: {projdir}  OK...')

        # move all save files into existing <tempdir>/{outfile_prefix}Solo.out dir. 
        for ext in MOVES:
            try:
                srcfile = f'{outfile_prefix}{ext}'
                destdir = f'{outfile_prefix}Solo.out'
                self.log.debug(f'moving {srcfile} -> {destdir} ...')
                shutil.move(srcfile, destdir)
            except FileNotFoundError:
                pass
        
        # move Solo.out dir to <cachedir> 
        outdir = f'{outfile_prefix}Solo.out'
        base = os.path.basename(outdir)
        dirname = os.path.dirname(outdir)
        self.log.debug(f'Got base of {base} dirname {dirname}')
        destdir = f'{projdir}/{base}'
        self.log.debug(f'destination directory name is {destdir}')
        newdir = shutil.copytree(outdir, destdir, dirs_exist_ok=True)
        # dst, symlinks, ignore, copy_function, ignore_dangling_symlinks, dirs_exist_ok)
        self.log.info(f'<tempdir> output copied to {newdir}')
        
        # clean tempdir. 
        self.log.debug(f'cleaning temp directory, removing {outdir}')
        shutil.rmtree(outdir)
        self.log.info(f'cleared temp dir of {outdir}')

    def _remove_fastqs(self, runlist):
        """
        Takes list of run_ids and removes <tempdir>/<run_id>*.fastq
        """
        self.log.debug(f'clearing tempdir of fastqs from {len(runlist)} runs...')
        for run_id in runlist:
            for fqfile in glob.glob(f'{self.tempdir}/{run_id}*.fastq'):
                os.remove(fqfile)
                self.log.debug(f'removed tempfile: {fqfile}')
         
        

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
