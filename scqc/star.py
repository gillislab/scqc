#!/usr/bin/env python
#
#  Module to deal with running STAR/StarSolo
#

import argparse
import glob
import importlib
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
from scqc.common import *


# inputs should be runs  identified as 'some10x'
# .sra files should already downloaded.
# run_id and species needs to be passed in from the dataframe
# can pass star parameters from config


# TODO queue 

class UnsupportedTechnologyException(Exception):
    """
    Thrown when run technology is neither 10x nor SmartSeq
    """

class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """
    
class NoTechnologyProjectException(Exception):
    """
    Thrown when a project has no runs with known technology.
    
    """


class Analyze(object):
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
        self.nocleanup = self.config.getboolean('star','nocleanup')
        self.force = self.config.getboolean('star','force')
        self.unzip = self.config.get('star','unzip')
        backstr = [ x.strip() for x in self.config.get('star','backends').split(',') ]    
        self.backends = {}
        for be in backstr:
            self.backends[be] = importlib.import_module(f'scqc.{be}')
        self.log.debug(f'STAR AlignReads initted. backends = {self.backends}')
       

    def execute(self, proj_id):
        # get relevant metadata
        proj_idf = self._get_meta_data(proj_id)
        self.log.debug(f'got impute for proj:\n{proj_idf}')
        proj_idf = self._known_tech(proj_idf)
        runlist  = list( proj_idf[ proj_idf.proj_id==proj_id ].run_id.unique() )
        self.log.debug(f'initializing STAR alignment for {proj_id} {len(runlist)} runs with known tech.')
        # should contain proj_id if category applies
        done = None
        part = None
        seen = proj_id

        # Overall flags. 
        somedone = False
        partial = False
        somefailed = False        
        # bring in all fastqs possible to <tempdir>
        # get first data_source value for project id. (assuming all are same/correct)
    
        try:
            backstr = proj_idf[ proj_idf.proj_id == proj_id].data_source.values[0]
            self.log.debug(f'got backend {backstr} for project {proj_id} ')
            self.backends[backstr].stage_in(self.config, self.cachedir, self.tempdir, runlist, self.force)
            
            # Overall flags. 
            somedone = False
            partial = False
            somefailed = False
                
            # split by technology and parse independently based on tech
            
            for tech, df in proj_idf.groupby(by = "tech_version") :
                self.log.debug(f'handling df=\n{df} with tech {tech}')
                try:
                    if tech =="smartseq":
                        # somedone, somefailed
                        (some, part) = self._handle_smartseq(proj_id, df)
                        if some:
                            # for smartseq, somedone means successful
                            somedone = True
                        else :
                            somefailed = True
                        self.log.debug(f'{proj_id} smartseq somedone={somedone} somefailed={somefailed}')
                            
                    elif tech.startswith('10xv'):  
                        (some, failed) = self._handle_10x(proj_id, tech, df)
                        if some:
                            # for 10x, some means some, possibly all
                            somedone = True
                        if failed:
                            somefailed = True
                        if some and failed:
                            partial = True
                        self.log.debug(f'{proj_id} {tech} somedone={somedone} somefailed={somefailed} partial={partial}')
    
                except Exception as ex:
                    self.log.error(f'fatal problem with NCBI proj_id {proj_id}')
                    self.log.error(traceback.format_exc(None))
                    
                finally:
                    # finally, clean up fastq files 
                    if not self.nocleanup:
                        self._cleantemp(proj_id, runlist)
                        self._remove_fastqs(proj_id, runlist)
                    else:
                        self.log.info(f'nocleanup is true. leaving temp files.')
        except Exception as ex:
            logging.error(f'something went wrong with {proj_id}')
            self.log.error(traceback.format_exc(None))
            
        if somedone and not somefailed:
            done = proj_id
            part = None
        elif somedone and somefailed:
            done = None
            part = proj_id
        else:
            done = None
            part = None
        
        self.log.debug(f'{proj_id} final: done={done} part={part} seen={seen}')
        return (done, part, seen)


    def _known_tech(self, df):
        """
        Filters df by known tech. 
        """
        KNOWN = ['smartseq','10xv3', '10xv2','10xv1']
        self.log.debug(f'filter by known tech. inlength={len(df)}')
        retdf = df[ df.tech_version.isin(KNOWN) ]
        self.log.debug(f'filter by known tech. outlength={len(retdf)}')        
        return retdf


    def _handle_smartseq(self, proj_id, df):
        # build the manifest
        df.reset_index(inplace=True, drop=True)
        (manipath, manifest) = self._make_manifest(proj_id, df) 
        # run star
        self.log.debug(f'Starting smartseq alignment for {proj_id}')
        try:
            outfile_prefix = self._run_star_smartseq(proj_id, manipath)
            self.log.debug(f'Got outfile_prefix={outfile_prefix} for {proj_id}')
            self._stage_out(proj_id, outfile_prefix)
            # (somedone, somefailed) 
            return (True, False)

        except Exception as ex:
            self.log.warning(f'smartseq run failed. ')
            self.log.error(traceback.format_exc(None))
            return (False, False)


    def _handle_10x(self, proj_id, tech, df):
        """
        For 10x, detect lanes. check read lengths, and trigger star run...
        Each invocation has unique tech (10xv2 or 10xv3)
        #for row in range(df.shape[0]):
        """
        partial = False
        somedone = False
        somefailed = False
        
        df.reset_index(inplace=True, drop=True)
        self.log.debug(f'df=\n{df}\nshape={df.shape}')
        
        df['prefix'] = df.apply(apply_striplane, axis=1)             
        for prefix, tdf in df.groupby(by = "prefix") :
            self.log.debug(f'starting {tech} processing for {prefix}')                
            try:
                cdnas = []
                barcodes = []
                gzipped = False
                for row in tdf.iterrows():
                    self.log.debug(f'row is {row}')
                    read1 = f"{self.tempdir}/{row[1]['read1']}" # biological cDNA
                    read2 = f"{self.tempdir}/{row[1]['read2']}" # technical CBarcode + UMI
                    cdnas.append(read1)
                    barcodes.append(read2)
                    if read1.endswith('.gz'):
                        gzipped = True
                # are there multiple lanes (i.e. multiple rows?)
                if len(tdf) > 1:
                    self.log.info(f'handling multiple lanes for prefix={prefix}')
                    cdnasarg = ','.join(cdnas)
                    barcodesarg = ','.join(barcodes)
                else:
                    self.log.info(f'handling single lane')
                    cdnasarg = cdnas.pop() 
                    barcodesarg = barcodes.pop()
                
                # check parity and run STAR
                # Short circuit if valid output exists.                 
                if len(cdnas)  == len(barcodes):
                    if not self._check_run_done(prefix, tech):
                        outfile_prefix = self._run_star_10x(prefix, tech, cdnasarg, barcodesarg, gzipped)
                        self.log.debug(f'Got outfile_prefix={outfile_prefix} for {proj_id} and {prefix}')
                        self._stage_out(proj_id, outfile_prefix)
                        self.log.debug(f'Stageout complete for {outfile_prefix}Solo.out')
                        somedone = True
                    else:
                        somedone = True
                else:
                    self.log.warning(f'mismatched cdnas/barcodes: {cdnas} <-> {barcodes} ')
            
            except Exception as ex:
                self.log.warning(f'Problem with run_id {prefix}.')
                self.log.error(traceback.format_exc(None))
                somefailed = True
            
            finally:
                if not self.nocleanup:
                    self._cleanrun(prefix)
                
                        
        self.log.info(f'completed handling for project {proj_id} tech={tech}')
        return( somedone, somefailed )


    def _check_run_done(self):
        return False

    
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
        idf = load_df(f'{self.metadir}/impute.tsv')
        # filter to include only requested project id
        proj_idf = idf[ idf.proj_id == proj_id]
        #impute = impute.loc[ impute.proj_id==proj_id ,:]
        # filter to include only requested species and only keep run ids
        #impute = impute.loc[ impute.taxon == int(spec_to_taxon(self.species)) ,:].reset_index(drop=True)
        return(proj_idf)


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


    def _run_star_smartseq(self, proj_id, manipath, gzipped= False):
        
        df = pd.read_csv(manipath, sep='\t')
        df.columns = ['read1','read2','run_id']
        fn = df['read1'][0]
        if fn.endswith('gz'):
            gzipped = True
        self.log.debug(f'auto-detected gzipped files in manifest.')
        
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

        if gzipped : 
            cmd +=  ['--readFilesCommand','zcat']
            #cmd +=  ['--readFilesCommand','gunzip -c']            
            
        run_command(cmd)
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
    def _run_star_10x(self, run_id, tech, bio_readpath, tech_readpath, gzipped = False):
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
        if gzipped :
            cmd += ['--readFilesCommand',f'{self.unzip}']
        try:
            run_command(cmd)
        except NonZeroReturnException as nzre:
            self.log.error(traceback.format_exc(None))
        
        finally:
            if not self.nocleanup:              
                self._cleanrun(run_id)
        
        return(outfile_prefix)


    def _stage_out(self, proj_id, outfile_prefix):
        """
        Stage out. 
        <tempdir>/{outfile_prefix}Solo.out
                  {outfile_prefix}Log.out
                  {outfile_prefix}Log.final.out
                  {outfile_prefix}SJ.out.tab
                  {outfile_prefix}manifest.tsv        
                  {outfile_prefix}Log.progress.out
                  
        for 10x
        f'{self.tempdir}/{run_id}_{tech}_Solo.out' -> SRR10285015_10xv2_Solo.out
            -> <cachedir>/proj_id/{run_id}_{tech}_Solo.out/
            
        for smartseq
        f'{self.tempdir}/{proj_id}_smartseq_Solo.out' -> SRP066963_smartseq_Solo.out/
            -> <cachedir>/proj_id/{proj_id}_{tech}_Solo.out/
        
        """
        
        MOVES = ['Log.out',
                 'Log.final.out',
                 'SJ.out.tab',
                 'manifest.tsv']
        DELETES = ['Log.progress.out']
        
        self.log.debug(f'called for project {proj_id} and outfile_prefix= {outfile_prefix}')
          
        # move all save files into existing <tempdir>/{outfile_prefix}Solo.out dir. 
        for ext in MOVES:
            try:
                srcfile = f'{outfile_prefix}{ext}'
                destdir = f'{outfile_prefix}Solo.out/'
                self.log.debug(f'moving {srcfile} -> {destdir} ...')
                shutil.move(srcfile, destdir)
            except FileNotFoundError:
                pass

        for ext in DELETES:
            try:
                srcfile = f'{outfile_prefix}{ext}'
                self.log.debug(f'removing {srcfile}...')
                os.remove(srcfile)
            except FileNotFoundError:
                pass

        projdir = f'{self.cachedir}/{proj_id}/'
        try:    # make project specific solo out directories. 
            os.makedirs(projdir)
            self.log.debug(f'created cache project dir {projdir}')
        except FileExistsError:
            self.log.warning(f'cache project dir already exists: {projdir}  OK...')

        
        # move Solo.out dir to <cachedir> 
        outdir = f'{outfile_prefix}Solo.out'
        base = os.path.basename(outdir)
        dirname = os.path.dirname(outdir)
        self.log.debug(f'Got base of {base} dirname {dirname}')
        destdir = f'{projdir}{base}'
        self.log.debug(f'outdir is {outdir} destdir is {destdir}')
        newdir = shutil.copytree(outdir, destdir, dirs_exist_ok=True)
        # dst, symlinks, ignore, copy_function, ignore_dangling_symlinks, dirs_exist_ok)
        self.log.debug(f'copy done. changing permissions for {newdir}')
        chmod_recurse('{newdir}')
        self.log.info(f'<tempdir> output copied to {newdir} and permissions adjusted.')
        
        # clean tempdir. 
        if not self.nocleanup:
            self.log.debug(f'cleaning temp directory, removing {outdir}')
            shutil.rmtree(outdir)
            self.log.info(f'cleared temp dir of {outdir}')
        else:
            self.log.info(f'nocleanup true. leaving files.')


    def _cleanrun(self, run_id):
        """
        Temp cleanup for each individual run of STAR. Mostly important when returns non-zero 
        due to error. 
        E.g.
        <run_id>_10xv3__STARtmp
        <run_id>_L001
        <run_id>_L002
        
        """
        self.log.debug(f'beginning run cleanup for runid/prefix: {run_id}')
        filedirlist =  glob.glob(f'{self.tempdir}/{run_id}_*10xv*')
        self.log.debug(f'files/dirs to delete: {filedirlist}')
        remove_pathlist(filedirlist)

    
    def _cleantemp(self, proj_id, runlist):
        """
        Belt-and-suspenders cleanup after whole project run. 
        """
        filedirlist =  glob.glob(f'{self.tempdir}/{proj_id}_smartseq_*')
        remove_pathlist(filedirlist)
        for rid in runlist:
            filelist = glob.glob(f'{self.tempdir}/{rid}_*10xv*')
            remove_pathlist(filelist)
        

    def _remove_fastqs(self, proj_id, runlist):
        """
        Takes list of run_ids and removes <tempdir>/<run_id>*.fastq
        """
        self.log.debug(f'clearing tempdir of fastqs from {len(runlist)} runs...')
        for run_id in runlist:
            for fqfile in glob.glob(f'{self.tempdir}/{run_id}*.fastq'):
                os.remove(fqfile)
                self.log.debug(f'removed tempfile: {fqfile}')

def apply_striplane(row):
    """ 
    """
    return re.sub('_L00[0-9]','',row['run_id'])
        

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
               "--genomeSAindexNbases", "12",  # for low memory
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
