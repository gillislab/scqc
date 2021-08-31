#!/usr/bin/env python
#
#  Module to deal with interactions with SRA and parsing SRA metadata.
#
# http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=
#
# Could use  SRR14584407 SRR14584408 in example..
# John Lee's working examples - SRP126648 for Smartseq and SRP243446 for 10xv2

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

from pathlib import Path
from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty
from requests.exceptions import ChunkedEncodingError

import xml.etree.ElementTree as et
import datetime as dt
import pandas as pd
import numpy as np

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *

# Translate between Python and SRAToolkit log levels for wrapped commands.
#  fatal|sys|int|err|warn|info|debug
LOGLEVELS = {
    10: 'debug',
    20: 'info',
    30: 'warn',
    40: 'err',
    50: 'fatal',
}

PROJ_COLUMNS = ['proj_id', 'ext_ids', 'title', 'abstract', 'submission_id']

SAMP_COLUMNS = ['samp_id', 'ext_ids',  'taxon',
                'sciname', 'title', 'attributes', 'proj_id', 'submission_id']

EXP_COLUMNS = ['exp_id', 'ext_ids',  'strategy',
               'source', 'lcp', 'samp_id', 'proj_id', 'submission_id']

RUN_COLUMNS = ['run_id', 'ext_ids', 'tot_spots', 'tot_bases', 'run_size', 'publish_date',
               'taxon', 'organism', 'nreads',  'basecounts', 'file_url','file_size','exp_id', 'samp_id', 'proj_id', 
               'submission_id' ]

IMPUTE_COLUMNS = ['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon','batch']


TECH_RES = {
    '10x'   : re.compile("10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'|Kit v1|PN-120233|10X V1", re.IGNORECASE),
    #'10xv1' : re.compile("", re.IGNORECASE),
    #'10xv2' : re.compile("v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2", re.IGNORECASE),
    #'10xv3' : re.compile("v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078", re.IGNORECASE),
    'smartseq' : re.compile("Smart-Seq|SmartSeq|Picelli|SMART Seq", re.IGNORECASE),
    'smarter' : re.compile("SMARTer", re.IGNORECASE),
    'dropseq' : re.compile("Cell 161, 1202-1214|Macosko|dropseq|drop-seq", re.IGNORECASE),
    'celseq'  : re.compile("CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq", re.IGNORECASE),
    'sortseq' : re.compile("Sort-seq|Sortseq|Sort seq", re.IGNORECASE),
    'seqwell' : re.compile("Seq-Well|seqwell", re.IGNORECASE),
    'biorad'  : re.compile("Bio-Rad|ddSeq", re.IGNORECASE),
    'indrops' : re.compile("inDrop|Klein|Zilionis", re.IGNORECASE),
    'marsseq2': re.compile("MARS-seq|MARSseq|Jaitin et al|jaitin", re.IGNORECASE),
    'tang'    : re.compile("Tang", re.IGNORECASE),
    # 'TruSeq':re.compile("TruSeq", re.IGNORECASE),
    'splitseq': re.compile("SPLiT-seq", re.IGNORECASE),
    'microwellseq': re.compile("Microwell-seq", re.IGNORECASE)
}

#deprecated
keywords = {
    "is10x": "10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'",
    "v3": "v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078",
    "v2": "v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2",
    "v1": "Kit v1|PN-120233|10X V1",
    "ss": "Smart-Seq|SmartSeq|Picelli|SMART Seq",
    "smarter": "SMARTer",
    "dropseq": "Cell 161, 1202-1214|Macosko|dropseq|drop-seq",
    "celseq": "CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq",
    "sortseq": "Sort-seq|Sortseq|Sort seq",
    "seqwell": "Seq-Well|seqwell",
    "biorad": "Bio-Rad|ddSeq",
    "indrops": "inDrop|Klein|Zilionis",
    "marsseq2": "MARS-seq|MARSseq|Jaitin et al|jaitin",
    "tang": "Tang",
    # "TruSeq":"TruSeq",
    "splitseq": "SPLiT-seq",
    "microwellseq": "Microwell-seq"
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


class RunUnavailableException(Exception):
    """ Thrown when Run in a Runset is unavailable.  """


class SampleUnavailableException(Exception):
    """ Thrown when Sample in a Runset is unavailable.  """

class MissingReadsException(Exception):
    """ Thrown when Run is lacking one or more reads.  """


def setup(config):
    '''
    Builds directories in config file 
    Download the appropriate supplement data.
    Only needs to be done once.
    '''

    log = logging.getLogger('setup')
    config = config
    # directories
    metadir = os.path.expanduser(config.get('setup', 'metadir'))
    cachedir = os.path.expanduser(config.get('setup', 'cachedir'))
    tempdir = os.path.expanduser(config.get('setup', 'tempdir'))
    resourcedir = os.path.expanduser(config.get('setup', 'resourcedir'))
    outputdir = os.path.expanduser(config.get('setup', 'outputdir'))
    figuredir =os.path.expanduser(config.get('setup', 'figuredir'))

    dirs_to_make = [metadir, 
                    f'{cachedir}/sra', 
                    tempdir, 
                    resourcedir, 
                    outputdir, 
                    figuredir]

    for direc in dirs_to_make:
        try:
            os.makedirs(direc)
        except FileExistsError:
            pass


# To do: batch for large queries.
# Ideally, we'd query once for all current data, (one large query)
# then periodically query every few days (many smaller queries)
# To do: Batch the requests to efetch
# looping through all IDs in efetch can be slow (?)
# requires more requests to the server and fails more often - my tests:  fails ~25% of the time
# Note: repeated iterations of Query.execute() will include the failed UIDs
class Query(object):
    """

    Run info for sample and projects?
    wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRS049712'
    wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP290125'

    """

    def __init__(self, config):
        self.log = logging.getLogger('sra')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('query', 'metadir'))
        self.cachedir = os.path.expanduser(
            self.config.get('query', 'cachedir'))
        self.species = self.config.get('sra','species')
        self.tissue = self.config.get('sra','tissue')
        self.sra_esearch = self.config.get('sra', 'sra_esearch')
        self.sra_efetch = self.config.get('sra', 'sra_efetch')
        self.search_term = self.config.get('sra', f'{self.species}_search_term')
        self.query_max = self.config.get('sra', 'query_max')
        self.expid_file = os.path.expanduser(self.config.get('sra', 'expid_file'))
        self.expids = self.build_expidset()
        self.query_sleep = float(self.config.get('sra', 'query_sleep'))
        self.xid_batchsize = int(self.config.get('sra', 'xid_batchsize'))


    def build_expidset(self):
        self.log.debug(f'expid_file is {self.expid_file}')
        expidset = set(readlist(self.expid_file))
        self.log.debug(f'expidset len={len(expidset)}')
        return expidset

    def execute(self, proj_id):
        """
        For proj_id:
            For all experiments that *are* in exppid file (omitting experiments in project 
            that didn't match a priori query). 
                Perform full query, get runids, fetch for each id, parse XML response. 
                Put project and run info in project_metadata.tsv and project_runs.tsv
                Put completed project ids into query-donelist.txt

        """
        self.log.info(f'handling proj_id {proj_id}')
        try:
            pdf = query_project_metadata(proj_id)
            explist = list(pdf.Experiment)
            self.log.info(
                f'proj_id {proj_id} has {len(explist)} experiments.')
            proj_rows = []
            samp_rows = []
            exp_rows = []
            run_rows = []
            
            # XXX batched this - JL
            # intersect explist (project specific) and global expids 
            explist = list(set(explist) & self.expids)
            nexps = len(explist)
            nbatches = int(np.ceil(nexps / self.xid_batchsize ))
            for i in range(nbatches): 
                exps = explist[ i*self.xid_batchsize: (i+1)*self.xid_batchsize  ]
                exps = ",".join(exps)
                # if exp in self.expids:
                exd = self.query_experiment_package_set(exps)
                (projrows, samprows, exprows,
                    runs) = self.parse_experiment_package_set(exd)
                proj_rows = itertools.chain(proj_rows, projrows)
                samp_rows = itertools.chain(samp_rows, samprows)
                exp_rows = itertools.chain(exp_rows, exprows)
                run_rows = itertools.chain(run_rows, runs)
                # else:
                #     self.log.debug(f'expid {exp} not in set of search expids.  ')
            proj_rows = list(proj_rows)
            samp_rows = list(samp_rows)
            exp_rows = list(exp_rows)
            run_rows = list(run_rows)

            self.log.debug(f'final proj_rows: {proj_rows}')
            self.log.debug(f'final samp_rows: {samp_rows}')
            self.log.debug(f'final exp_rows: {exp_rows}')
            self.log.debug(f'final run_rows: {run_rows}')
            
            # make dataframes
            pdf = pd.DataFrame(proj_rows, columns=PROJ_COLUMNS)
            sdf = pd.DataFrame(samp_rows, columns=SAMP_COLUMNS)
            edf = pd.DataFrame(exp_rows, columns=EXP_COLUMNS)
            rdf = pd.DataFrame(run_rows, columns=RUN_COLUMNS)
            
            # merge dataframes to files. 
            merge_write_df(pdf, f'{self.metadir}/projects.tsv')            
            merge_write_df(sdf, f'{self.metadir}/samples.tsv')
            merge_write_df(edf, f'{self.metadir}/experiments.tsv')
            merge_write_df(rdf, f'{self.metadir}/runs.tsv')

            self.log.info(f'successfully processed project {proj_id}')
            # return proj_id only if it has completed successfully.
            return (proj_id, proj_id, proj_id)

        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            logging.error(traceback.format_exc(None))
            return(None, None, proj_id)


    def query_experiment_package_set(self, xid):
        """
        Query XML data for this experiment ID. 

        """
        xmldata = None
        try:
            url = f"{self.sra_efetch}&id={xid}"
            self.log.debug(f"fetch url={url}")

            # XXX could run infinitely???
            while True:
                r = requests.post(url)
                if r.status_code == 200:
                    xmldata = r.content.decode()
                    self.log.debug(f'good HTTP response for {xid}')
                    break
                else:
                    self.log.warn(
                        f'bad HTTP response for id {xid}. retry in 10s')
                    time.sleep(10)

        except Exception as ex:
            self.log.error(f'problem with NCBI id {xid}')
            logging.error(traceback.format_exc(None))

        finally:
            self.log.debug(
                f"sleeping {self.query_sleep} secs between fetch calls...")
            time.sleep(self.query_sleep)
        return xmldata


    def parse_experiment_package_set(self, xmlstr):
        """
        package sets should have one package per uid pulled via efetch, e.g.

        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=12277089,12277091
        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=13333495 

        """
        root = et.fromstring(xmlstr)
        self.log.debug(f"root={root}")
        proj_rows = []
        samp_rows = []
        exp_rows = []
        run_rows = []

        n_processed = 0
        for exp in root.iter("EXPERIMENT_PACKAGE"):
            (projrow, samprow, exprow, newruns) = self.parse_experiment_package(exp)
            proj_rows.append(projrow)
            samp_rows.append(samprow)
            exp_rows.append(exprow)
            # run_rows.append(newruns)
            run_rows = itertools.chain(run_rows, newruns)
            n_processed += 1
        self.log.debug(f"processed {n_processed} experiment package(s).")
        run_rows = list(run_rows)
        self.log.debug(
            f'returning\n    proj_rows: {proj_rows}\n    exp_rows: {exp_rows} \n    run_rows: {run_rows}')
        return (proj_rows, samp_rows, exp_rows, run_rows)

    def parse_experiment_package(self, root):
        """
        NCBI provides no XSD, so we shouldn't rely on order

        """
        self.log.debug('parsing experiment package...')
        exp = root.find('EXPERIMENT')
        sub = root.find('SUBMISSION')
        proj = root.find('STUDY')
        samp = root.find('SAMPLE')
        runs = root.find('RUN_SET')

        # get experiment properties.

        # get submission properties
        sra_id = sub.get('accession')

        # get study/project properties title, abstract
        projrow = self.parse_proj(proj)
        projrow.append(sra_id)
        proj_id = projrow[0]

        # get sample properties - append project and sra ids
        samprow = self.parse_sample(samp)
        samprow.append(proj_id)
        samprow.append(sra_id)
        samp_id = samprow[0]

        # get experiment properties
        exprow = self.parse_exp(exp)
        exprow.append(sra_id)
        exp_id = exprow[0]

        # get run properties - list of lists
        runrows = self.parse_run_set(runs, proj_id, sra_id)

        self.log.debug(
            f'exprow: proj_id={proj_id} exp_id={exp_id} sra_id={sra_id} samp_id={samp_id}')
        # projrow = [proj_id, title, pubdate, abstract]
        # exprow = [proj_id, exp_id, sra_id, gsm, gse, lcp, sample_attributes]
        self.log.debug(
            f'\n  projrow: {projrow}\n   exprow: {exprow} \n  runrows: {runrows}')
        return(projrow, samprow, exprow, runrows)

    def parse_run_set(self, runs, proj_id, sra_id):
        """

        """
        runrows = []
        for run in runs.findall('RUN'):
            runrow = self.parse_run(run)
            runrow.append(proj_id)
            runrow.append(sra_id)
            runrows.append(runrow)
        return runrows

    # need to append project id
    def parse_run(self, run):

        run_ext_ids = {}
        ids = run.find('IDENTIFIERS')
        run_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            run_ext_ids[tag] = val
        run_ext_ids = str(run_ext_ids)
        avail_status = run.get('unavailable')
        if avail_status == 'true':
            raise RunUnavailableException(f'run data unavailable for {run_id}')

        total_spots = run.get('total_spots')
        total_bases = run.get('total_bases')
        run_size = run.get('size')
        pdate = run.get('published')

        expid = run.find('EXPERIMENT_REF').get('accession')
        pool = run.find('Pool')
        sampleid = pool.find('Member').get('accession')
        taxon = pool.find('Member').get('tax_id')
        organism = pool.find('Member').get('organism')

        nreads = run.find('Statistics').get('nreads')

        srafiles = run.find('SRAFiles')

        (file_url, file_size) = self.parse_srafiles(srafiles, run_id)

        bases = run.find('Bases')
        basecounts = {}
        for base in bases:
            tag = base.get('value')
            val = base.get('count')
            basecounts[tag] = val

        basecounts = str(basecounts)

        runrow = [run_id, run_ext_ids, total_spots, total_bases, run_size, pdate,
                  taxon, organism, nreads,  basecounts, file_url, file_size , expid, sampleid ]

        return runrow

    def parse_srafiles(self, srafiles, runid):
        """
        Get file info. Choose Amazon URL if available. 
        
        """
        file_url = None
        file_size = None
        for srafile in srafiles.findall('SRAFile'):
            if srafile.get('filename') == runid:
                file_size = srafile.get('size')
                file_url = srafile.get('url')
                for altern in srafile.findall('Alternatives'):
                    url = altern.get('url')
                    if 'amazonaws.com' in url:
                        self.log.debug(f'found AWS alternate: {url}  Using...')
                        file_url = url
        self.log.info(f'got info for {runid}: {file_size} {file_url}')
        return (file_url, file_size)

    def parse_exp(self, exp):

        exp_ext_ids = {}
        ids = exp.find('IDENTIFIERS')
        exp_id = ids.find('PRIMARY_ID').text
        self.log.debug(f'parsing exp_id: {exp_id}')
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            exp_ext_ids[tag] = val
        exp_ext_ids = str(exp_ext_ids)
        projid = exp.find('STUDY_REF').get('accession')
        des = exp.find('DESIGN')
        sampid = des.find('SAMPLE_DESCRIPTOR').get('accession')
        ldes = des.find('LIBRARY_DESCRIPTOR')

        lcp = ""
        strat = ""
        source = ""
        try:
            lcp = ldes.find('LIBRARY_CONSTRUCTION_PROTOCOL').text
            lcp = lcp.strip()
            strat = ldes.find('LIBRARY_STRATEGY').text
            source = ldes.find('LIBRARY_SOURCE').text
        except AttributeError as ae:
            self.log.warn(f'attribute error parsing LIBRARY_DESCRIPTOR children.')

        exprow = [exp_id, exp_ext_ids,  strat, source, lcp,  sampid, projid]
        return exprow

    # need to append project id - done in execute
    def parse_sample(self, samp):
        samp_ext_ids = {}
        ids = samp.find('IDENTIFIERS')
        samp_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            samp_ext_ids[tag] = val
        samp_ext_ids = str(samp_ext_ids)
        
        try:
            samptitle = samp.find('TITLE').text
        except:
            samptitle = samp_id

        sample_attributes = {}
        for elem in samp.find('SAMPLE_ATTRIBUTES').findall('SAMPLE_ATTRIBUTE'):
            tag = elem.find('TAG').text
            val = elem.find('VALUE').text
            sample_attributes[tag] = val
        sample_attributes = str(sample_attributes)

        taxid = samp.find('SAMPLE_NAME').find('TAXON_ID').text
        sciname = samp.find('SAMPLE_NAME').find('SCIENTIFIC_NAME').text
        samprow = [samp_id, samp_ext_ids,  taxid,
                   sciname, samptitle, sample_attributes]

        return samprow

    def parse_proj(self, proj):

        proj_ext_ids = {}
        ids = proj.find('IDENTIFIERS')
        proj_id = ids.find('PRIMARY_ID').text
        for elem in ids.findall('EXTERNAL_ID'):
            tag = elem.get('namespace')
            val = elem.text
            proj_ext_ids[tag] = val
        proj_ext_ids = str(proj_ext_ids)  # convert to strings to store in df
        d_elem = proj.find('DESCRIPTOR')
        title = d_elem.find('STUDY_TITLE').text
        abstract = d_elem.find('STUDY_ABSTRACT').text

        projrow = [proj_id, proj_ext_ids, title, abstract]
        return projrow

    def query_runs_for_project(self, project):
        """     
        wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP290125'
            Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
            SRR12951718,2020-11-05 13:15:13,2020-11-05 12:55:40,128538925,12596814650,0,98,3781,GCA_000001635.4,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR12951718/SRR12951718.1,SRX9404801,,RNA-Seq,cDNA,TRANSCRIPTOMIC,SINGLE,0,0,ILLUMINA,NextSeq 550,SRP290125,PRJNA673364,3,673364,SRS7622575,SAMN16604770,simple,10090,Mus musculus,GSM4873966,,,,,,,no,,,,,GEO,SRA1151197,,public,1790FB1FF1C3B1A1D0E1958BE6859830,86BE0A6ECACE94B3BC53AFA6FE2A258F
            SRR12951719,2020-11-05 12:52:07,2020-11-05 12:39:41,137453038,13470397724,0,98,4581,GCA_000001635.4,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR12951719/SRR12951719.1,SRX9404802,,RNA-Seq,cDNA,TRANSCRIPTOMIC,SINGLE,0,0,ILLUMINA,NextSeq 550,SRP290125,PRJNA673364,3,673364,SRS7622574,SAMN16604769,simple,10090,Mus musculus,GSM4873967,,,,,,,no,,,,,GEO,SRA1151197,,public,8C05F51EF726335BE2C41DB9A6323EC4,2C810A53A72F30ECAC8AAD71B4E55CAB
            SRR12951720,2020-11-05 12:26:06,2020-11-05 12:13:08,100175297,9817179106,0,98,2840,GCA_000001635.4,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR12951720/SRR12951720.1,SRX9404803,,RNA-Seq,cDNA,TRANSCRIPTOMIC,SINGLE,0,0,ILLUMINA,NextSeq 550,SRP290125,PRJNA673364,3,673364,SRS7622576,SAMN16604768,simple,10090,Mus musculus,GSM4873968,,,,,,,no,,,,,GEO,SRA1151197,,public,D02C61154EB828AB07968FF1BFE52485,D80E2EB7A3A50FCE5062F584FD58FD8F


        """
        pass

    def _split_df_by_project(self, df):
        self.metadir
        for srp, srp_df in df.groupby('project', as_index=False):

            for tech, srp_tech_df in srp_df.groupby('method', as_index=False):
                outfile = f'{self.metadir}/{srp}_metadata.tsv'
                srp_tech_df.to_csv(outfile, sep="\t", mode="a",
                                   index=False, header=not os.path.exists(outfile))

        return




class Impute(object):
    """
    Imputes sequencing technology for all runs under a project. 

    """
    def __init__(self, config):
        self.log = logging.getLogger('impute')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('impute', 'metadir'))


    def execute(self, proj_id):
        """
        XXX if tech is not supported, will return an empty dataframe
        For proj_id:
            Impute technology where possible. 
            Put completed project ids into query-donelist.txt
            examples    - SRP114926 - contains both 10x and smartseq
                        - SRP122508 - contains just 10xv2. 192 runs, 10 exp, 10 samples
        """
        self.log.info(f'handling proj_id {proj_id}')
        try:
            # read in experiment file
            expfile = f'{self.metadir}/experiments.tsv'
            edf = load_df(expfile)
            self.log.debug(f'opened experiments DF OK...')
            edf = edf[edf.proj_id == proj_id].reset_index(drop=True) # rename pdf -> edf 
            self.log.debug(f'got project-specific df: \n{edf}')
            # impute technology  -  exp_id|tech
            idf = self.impute_tech_from_lcp(edf)    


            ss_manual = pd.read_csv(f'{self.metadir}/smartseq_projs.tsv' , header=None)
            ss_manual.columns = ['proj_id']
            if proj_id in ss_manual.proj_id.unique() :
                # read the vdb dump data
                vdf = load_df(f'{self.metadir}/vdb_dump.tsv')
                vdf = vdf[vdf.proj_id == proj_id].reset_index(drop=True) 
                # only take the consistent runs lengths 
                m = vdf.read_lengths.value_counts().index[0]
                print(vdf.read_lengths.value_counts())
                print(m)
                
                vdf = vdf.loc[vdf.read_lengths == m,:]
                #update idf
                idf.tech[idf.exp_id.isin(vdf.exp_id)] = 'smartseq'
            
            
            tx_manual = pd.read_csv(f'{self.metadir}/tenx_projs.tsv' , header=None)
            tx_manual.columns = ['proj_id']
            if proj_id in tx_manual.proj_id.unique() :
                
                vdf = load_df(f'{self.metadir}/vdb_dump.tsv')
                vdf = vdf[vdf.proj_id == proj_id].reset_index(drop=True) 
                # only take runs with read_length containing 24, 26, 28
                m = vdf.read_lengths.str.contains('24,|, 24|26,|, 26|28,|, 28')
                vdf = vdf.loc[m,:]
                
                print(vdf.read_lengths.value_counts())
                
                #update idf
                idf.tech[idf.exp_id.isin(vdf.exp_id)] = '10x'




            self.log.debug(f'got initial imputed tech df: \n{idf}')

            # match run to tech
            runfile = f'{self.metadir}/runs.tsv'
            rdf = load_df(runfile)
            rdf = rdf[rdf.proj_id == proj_id].reset_index(drop=True)
            # impute 10x version
                

            outdf = self.impute_10x_version(idf, rdf)
            self.log.debug(f'got imputed 10x version df: \n{outdf}')
            
            ssdf = self.parse_smartseq(idf, rdf)
            self.log.debug(f'parsed smartseq df: \n{ssdf}')
            # save to disk
            outdf=outdf.append(ssdf)

            # append the inferred batch from samples.tsv
            samplefile = f'{self.metadir}/samples.tsv'
            sdf = load_df(samplefile)
            sdf = sdf[sdf.proj_id == proj_id].reset_index(drop=True)

            #impute batch
            bdf = self.impute_batch(sdf, rdf)
            outdf = outdf.merge(bdf, on='run_id',how='inner')
            self.log.debug(f'batches inferred: \n{bdf}')
            # save to disk
            outdf = outdf[['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon','batch']]
            outdf.columns = IMPUTE_COLUMNS  # renames the columns from global 
            
            outdf = self._known_tech(outdf)
            if len(outdf) > 0:
                merge_write_df(outdf, f'{self.metadir}/impute.tsv')  
            else :
                self.log.warn(f'Unable to predict tech for:{proj_id} ')
                return (None, None, proj_id)
            
            self.log.info(f'completed imputation for proj_id {proj_id}')
            return (proj_id, proj_id, proj_id)          

        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            logging.error(traceback.format_exc(None))
            return (None, None, proj_id)


    def _known_tech(self, df):
        """
        Filters impute df by known tech. 
        """
        KNOWN = ['smartseq','10xv3', '10xv2','10xv1']
        self.log.debug(f'filter by known tech. inlength={len(df)}')
        retdf = df[ df.tech_version.isin(KNOWN) ]
        self.log.debug(f'filter by known tech. outlength={len(retdf)}')        
        return retdf
            
    # TODO deal with multiple tech finds
    def impute_tech_from_lcp(self, df):
        '''
        Take in experiment df for specific project.  
            Get unique library construction protocol (lcp) values. 
            For each value, loop over all keyword terms, identifying putative tech. 
            Fill in method column for all runs in temp DF.   
            Create impute DF
        
        '''
        logging.debug(f'got df: \n{df}')
        df.lcp= df.lcp.fillna('None').values

        ulcp = pd.DataFrame({"lcp": df.lcp.unique()})
        ulcp['tech'] ='unknown'
        ulcp['techcount'] =0
        # ulcp= ulcp.loc[1:,:].reset_index(drop=True) 
        
        # search for the keywords
        # doesn't play nice with NaN lcp - fill with a str

        for i in range(len(TECH_RES)):
            key = list(TECH_RES)[i]
            kw = TECH_RES[key]
            ulcp[key] = ulcp.lcp.str.contains(kw)
        
            ulcp.loc[ulcp[key],'tech'] = key
            ulcp.loc[ulcp[key],'techcount'] +=1
        
        ulcp.loc[ulcp['techcount'] > 1,'tech'] = 'multiple'
        dfout = ulcp[['lcp','tech']]
        self.log.debug(f'keyword hits: {ulcp.tech.values}')
        
        dfout = dfout.merge(df[['exp_id','lcp']], on="lcp")

        return dfout[['exp_id','tech']]



    # TODO error handling
    def impute_10x_version(self, idf, rdf):
        """
        For known 10x, get first part of fasta file and determine version.
        Only looks at the 10x portion if multiple techs used.
        Returns only 10x runs
        doesn't matter if multiple projects are included
        """

        # look at the first few lines of the fastq file.
        loglev = LOGLEVELS[self.log.getEffectiveLevel()]
        # print( f'\n{loglev}\n')
        
        # XXX require runs have rdf.nreads > 1. Otherwise, unable to impute
        #self.log.debug(f'idf={idf} rdf={rdf.nreads}')
        ind = rdf.nreads > 1
        self.log.debug(f'number of runs where nreads <= 1 {sum(ind)}') 
        rdf = rdf [ind]
        df = rdf.merge(idf, on = 'exp_id',how='left')

        # get all runs associated with the 10x inferred experiments
        runs = df.loc[ df.tech == '10x','run_id']
        if len(runs)  == 0:
            self.log.debug('no 10x runs imputable') 
           
            return pd.DataFrame(columns=['run_id' ,'tech_version','read1','read2','exp_id','proj_id', 'taxon'])

        allRows = []
        # TODO multithread?
        for run_id in runs :
            cmd = ['vdb-dump', 
                '--rows', '1',
                '--columns', 'READ_LEN',
                run_id]


            # TODO suppress messages from cmd
            cmdstr = " ".join(cmd)
            self.log.debug(f"vdb-dump command: {cmdstr} running...")
            cp = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            lengths = cp.stdout.read().decode("utf-8").strip().split(': ')[1].split(', ')

            self.log.debug(f"got {lengths}")
            ind = lengths.index(max(lengths))

            # which is the longest? use as cDNA read
            read_bio = f'{run_id}_{ind+1}.fastq'
            
            # 10xv2 is typically 98 bp
            # 10xv3 is typically 91 bp
            tech = "10x"
            for i in range(len(lengths)):
                if lengths[i] == '24':
                    tech = "10xv1"
                    ind2 = i
                    break
                elif lengths[i] == '26':
                    ind2 = i
                    tech = "10xv2"
                    break
                elif lengths[i] == '28':
                    ind2 = i
                    tech = "10xv3"
                    break


            if tech == "10x":
                read_tech = '-'
            else:
                read_tech = f'{run_id}_{ind2+1}.fastq'
            
            if read_tech == read_bio :
                read_tech = '-'
                tech = '10x'
            allRows.append([run_id, tech, read_bio, read_tech])

        outdf = pd.DataFrame(allRows,columns=['run_id','tech_version','read1','read2'] )
        tax = rdf[['run_id','samp_id','exp_id','proj_id','taxon']]

        outdf = outdf.merge(tax, on='run_id', how='inner') 
        
        # make sure to return in order!
        return outdf[['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon']]


    def parse_smartseq(self,idf,rdf):
        """
        For known smartseq, parse, build manifest
        doesn't matter if multiple projects are included
        """

        log = logging.getLogger('sra')
        # look at the first few lines of the fastq file.
        loglev = LOGLEVELS[log.getEffectiveLevel()]

        # get all runs associated with the smartseq inferred experiments
        df = rdf.merge(idf , on="exp_id", how = 'inner')    
        df = df[['run_id','taxon','nreads','exp_id','tech']]
        runs = df.loc[ df.tech == 'smartseq','run_id']
        
        outdf=pd.DataFrame({'run_id' : runs ,'tech_version':'smartseq'})

        # used to construct manifest and to be consistent with 10x section
        outdf.loc[df.nreads > 1 ,'read1']  = outdf['run_id']+ '_1.fastq'
        outdf.loc[df.nreads > 1 ,'read2']  = outdf['run_id']+ '_2.fastq'
        outdf.loc[df.nreads == 1 ,'read1'] = outdf['run_id']+ '.fastq'
        outdf.loc[df.nreads == 1 ,'read2'] = '-'
        tax = rdf[['run_id','exp_id','samp_id','proj_id','taxon']]

        outdf = outdf.merge(tax, on='run_id', how='inner') 
        
        # make sure to return in order!
        return outdf[['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon']]

    # sample attributes alone aren't enough to impute batch. use freq of ids
    def impute_batch(self, sdf, rdf):
        '''
        Uses the sample data to infer batch. If no batches are found, (i.e. everything 
        gets assigned batch 0), use cell/runs > as batch predictor during `gatherstats.py`

        This should only be run at the project level dataframes, but just in case,
        Splits by project id first and assigns a set batches to each 
        project id 
        '''
        
        cols = rdf.columns.tolist()
        cols.append('batch')
        new_rdf = pd.DataFrame(columns=cols )
        for proj , df in sdf.groupby('proj_id') :
            
            tm_rdf = rdf.loc[ rdf.proj_id == proj ,:] 
            
            # factorize outputs a tuple (index, attribute)
            samp2batch = pd.DataFrame({ 'samp_id' : df['samp_id'],
                            'samp' : pd.factorize(df.samp_id)[0] ,
                            'attr': pd.factorize(df['attributes'])[0],
                            'exp_id': pd.factorize(df['samp_id'])[0]})
            self.log.debug(f'sample to batch \n{samp2batch}')
            # batches should  contain at least two runs
            counts1 = samp2batch.iloc[:,1:].apply( lambda x: len(x.unique()) ,axis=0) > 1
            # batches shouldn't include all cells
            counts2 = samp2batch.iloc[:,1:].apply( lambda x: len(x.unique()) ,axis=0) < tm_rdf.shape[0]

            batch = list(counts1 &counts2)    # which columns denote possible batches
            # indicates that the 'samp_id' shouldn't be used for batch- only for merge
            batch.insert(0,False)   
            # batch not found at the sample level. Look at other levels.
            # samples are either unique per run or there is only one sample
            

            if sum(batch) == 0 :
                run2batch = pd.DataFrame({ 'run_id' : pd.factorize(tm_rdf['run_id'])[0] ,    
                            'exp_id': pd.factorize(tm_rdf['exp_id'])[0]})
                counts1 = run2batch.apply( lambda x: len(x.unique()) ,axis=0) > 1
                # batches shouldn't include all cells
                counts2 = run2batch.apply( lambda x: len(x.unique()) ,axis=0) < tm_rdf.shape[0]
                batch = counts1 &counts2    # which columns denote possible batches
                
                # batch is still False for everything. No batch inferred. 
                # If no batch inferred, need to align, get number of cells/run.
                # if ncells/run > 1 -> assign runs as the batch (e.g. 10x samples)

                if batch.sum() == 0:
                    samp2batch['batch'] = 0 # everything gets the same batch id
                else:
                    samp2batch['batch'] = run2batch.loc[:,batch].iloc[:,-1]
            else :
                    samp2batch['batch'] = samp2batch.loc[:,batch].iloc[:,-1]

            # merge these batches with the runs 
            tm_rdf=tm_rdf.merge(samp2batch[['samp_id','batch']], how ="left", on="samp_id")

            new_rdf = pd.concat([new_rdf, tm_rdf])

        return new_rdf[['run_id','batch']]

# do this in __main__ to allow for multiple threads 
class Download(object):
    """
    Handles all runid SRA file downloads for a project. 
    For the SRA case, that means getting SRR id for each project, and downloading them. 


    """

    def __init__(self, config): #, outlist
        self.log = logging.getLogger('sra')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('sra', 'metadir'))
        self.cachedir = os.path.expanduser(self.config.get('sra', 'cachedir'))
        self.dltool = os.path.expanduser(self.config.get('download', 'dltool'))
        self.max_rate = os.path.expanduser(self.config.get('download', 'max_rate'))
        rdf_file = f'{self.metadir}/runs.tsv'
        self.rdf = load_df(rdf_file).drop_duplicates()        
        self.log.info(f'Download initialized. ')


    def execute(self, proj_id):
        ddf = self.rdf[self.rdf.proj_id == proj_id ]
        rundict = dict(zip(ddf.run_id, ddf.file_url))
        runlist = list(ddf.run_id)
        self.log.debug(f'{len(runlist)} runs in project {proj_id}')
        totaldiskspace = ddf.file_size.sum() * 1e-9
        self.log.debug(f'Expected disk space for SRA files for {proj_id} is {totaldiskspace} GB')
        donelist = []        
        triedlist = []
        runlength = len(runlist)
       
        if not os.path.isdir(f'{self.cachedir}/sra'):
            self.log.debug(f'making sra cache subdir...')
            os.mkdir(f'{self.cachedir}/sra')
        
        try:    
            for i, (runid, srcurl) in enumerate(rundict.items()):  
                if self.dltool == 'wget':
                    self.log.info(f'tool is wget. handling file_urls')
                    self.log.debug(f'handling runid {runid} srcurl {srcurl}')
                    destpath = f'{self.cachedir}/sra/{runid}.sra'
                    rc = download_wget(srcurl, destpath, 
                                       finalname=None, overwrite=True, decompress=True, 
                                       rate=f'{self.max_rate}')
                    if str(rc) == '0' :
                        self.log.debug(f'runid {runid} [{i+1}/{runlength}] handled successfully.')
                        donelist.append(runid)
                    else:
                        self.log.debug(f'runid {runid} failed. rc={rc}')
                    triedlist.append(runid)
                    
            # determine if we succesfully completed all runs for project.
            diffset = set(donelist).difference(set(runlist))
            self.log.debug(f'diffset is {diffset}')
            if len(diffset) > 0:
                self.log.error(f'{len(donelist)} of {len(runlist)} downloaded for proj_id {proj_id} ')
                return (None, None, proj_id)
            else:
                self.log.info(f'download successful for proj_id {proj_id}')
                return (proj_id, proj_id, proj_id)    
        
        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            self.log.error(traceback.format_exc(None))
            return (None, None, proj_id)

# inputs are the runs completed by prefetch
# assumes path is cachedir/<run>.sra
class FasterqDump(object):
    '''
        Simple wrapper for NCBI fasterq-dump

        Usage: fasterq-dump [ options ] [ accessions(s)... ]
        Parameters:
            accessions(s)                list of accessions to process
        Options:
        -o|--outfile <path>              full path of outputfile (overrides usage
                                         of current directory and given accession)
        -O|--outdir <path>               path for outputfile (overrides usage of
                                         current directory, but uses given
                                         accession)
        -b|--bufsize <size>              size of file-buffer (dflt=1MB, takes
                                         number or number and unit where unit is
                                         one of (K|M|G) case-insensitive)
        -c|--curcache <size>             size of cursor-cache (dflt=10MB, takes
                                         number or number and unit where unit is
                                         one of (K|M|G) case-insensitive)
        -m|--mem <size>                  memory limit for sorting (dflt=100MB,
                                         takes number or number and unit where
                                         unit is one of (K|M|G) case-insensitive)
        -t|--temp <path>                 path to directory for temp. files
                                         (dflt=current dir.)
        -e|--threads <count>             how many threads to use (dflt=6)
        -s|--split-spot                  split spots into reads
        -S|--split-files                 write reads into different files
        -3|--split-3                     writes single reads into special file
        -f|--force                       force overwrite of existing file(s)
        --include-technical           explicitly include technical reads
        -v|--verbose                     Increase the verbosity of the program
                                         status messages. Use multiple times for
                                         more verbosity.

    '''

    def __init__(self, config, run_id): #, outlist
        self.log = logging.getLogger('sra')
        self.run_id = run_id

        self.log.debug(f'handling id {run_id}')
        self.config = config
        self.cachedir = os.path.expanduser(
            self.config.get('download', 'cachedir'))
        self.metadir = os.path.expanduser(
            self.config.get('download', 'metadir'))
        self.tempdir = os.path.expanduser(
            self.config.get('download', 'tempdir'))
        self.threads = self.config.get('sra', 'fq_nthreads')
        self.force = self.config.getboolean('download','force')
        self.nocleanup = self.config.getboolean('download','nocleanup')


    def execute(self):
        """
        fasterq-dump creates all files in a temp file before moving to final output. 
        if final output exists, files are done. 
                
        """
        if self._files_exist() and not self.force:
            self.log.info(f'output files already exist for {self.run_id}')
            return 0
        else:
            try:
                self.log.debug(f'running fasterq_dump {self.run_id}')
                self.run_fasterq_dump()
                return 0
                self.log.info(f'fasterq_dump run without exception {self.run_id}')
                
            except Exception as ex:
                self.log.warning(f'problem extracting for {self.run_id}')
     

    def run_fasterq_dump(self):
 
        self.log.debug(f'extracting id {self.run_id}')
        loglev = LOGLEVELS[self.log.getEffectiveLevel()]
        cmd = ['fasterq-dump', 
            '--split-files',
            '--include-technical',
            '--force', 
            '--threads', f'{self.threads}',
            '--outdir', f'{self.tempdir}/',
            '-t', f'{self.tempdir}/',
            '--log-level', f'{loglev}',
            f'{self.cachedir}/sra/{self.run_id}.sra']

        cmdstr = " ".join(cmd)
        logging.info(f"fasterq-dump command: {cmdstr} running...")
        start = dt.datetime.now()
        cp = subprocess.run(cmd, 
                        text=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
        end = dt.datetime.now()
        elapsed =  end - start
        logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
        if str(cp.returncode) == '0':
            logging.info(f'successfully extracted run {self.run_id}')
        else:
            logging.debug(f"got stderr: {cp.stderr}")
            logging.debug(f"got stdout: {cp.stdout}")
            logging.error(f' unknown non-zero return code for {self.run_id}')       
        
        if not self._files_exist():
            raise MissingReadsException(f'For {self.run_id}')                        
  

    def _handle_incomplete(self):
        """
        check /temp/fasterq.tmp.*/ for files for this runid. 
        delete if they exist
        return True if some were found.
        return False if not. 
        """
        found = False
        flist = glob.glob(f'{self.tempdir}/fasterq.tmp.*/{self.run_id}*')
        if len(flist) > 0:
            found = True
            fqtemp =   os.path.dirname(flist[0])
            self.log.debug(f'fasterq tempdir = {fqtemp}')
            self.log.debug(f'found incomplete files {flist}')
            if not self.nocleanup:
                self.log.debug(f'removing fasterq tempdir {fqtemp} for runid {self.run_id}')
                remove_pathlist([fqtemp])
        return found

    def _files_exist(self):
        """
        Determine if output files already exist.
        E.g. SRR13782545_1.fastq    
             SRR13782545_2.fastq               
             SRR13782545_3.fastq
             
              
        Return True if done
        False if not. 
        """
        found = False
        flist = glob.glob(f'{self.tempdir}/{self.run_id}_*.fastq')
        if len(flist) > 0:
            found = True
            self.log.debug(f'found output files for {self.run_id}: {flist}')
        else:
            self.log.debug(f'no output files for {self.run_id}')
        return found
        

def get_runs_for_project(config, proj_id):
    """
    
    """
    
    metadir = os.path.expanduser(config.get('sra', 'metadir'))
    filepath = f"{metadir}/runs.tsv"
    if os.path.isfile(filepath):
        pdf = pd.read_csv(filepath, sep='\t', index_col=0)
        return list(pdf[pdf.proj_id == proj_id].run_id)
    else:
        return []


def query_project_metadata(project_id):
    '''
    E.g. https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?db=sra&rettype=runinfo&save=efetch&term=SRP131661

    wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP131661'    

    '''
    log = logging.getLogger('sra')
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi"

    headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding": "gzip,deflate,sdch",
        "Accept-Language": "en-US,en;q=0.8",
        "Cache-Control": "no-cache",
        "Connection": "keep-alive",
        "DNT": "1",
        "Host": "trace.ncbi.nlm.nih.gov",
        "Origin": "http://trace.ncbi.nlm.nih.gov",
        "Pragma": "no-cache",
        "Referer": "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?db=sra",
        "User-Agent": "Mozilla/5.0 (iPhone; CPU iPhone OS 6_0 like Mac OS X) AppleWebKit/536.26 (KHTML, like Gecko) Version/6.0 Mobile/10A5376e Safari/8536.25"}

    payload = {
        "db": "sra",
        "rettype": "runinfo",
        "save": "efetch",
        "term": project_id}

    df = None
    log.debug('opening request...')
    r = requests.put(url, data=payload, headers=headers, stream=True)
    if r.status_code == 200:
        log.info('got good return. reading CSV to dataframe.')
        with io.BytesIO(r.content) as imf:
            df = pd.read_csv(imf)
        return df

    else:
        log.warning(f'bad HTTP return for proj: {project_id}')
        raise Exception(f'bad HTTP return for proj: {project_id}')



# to do: include drivers for 10x and ss alignments
if __name__ == "__main__":

    gitpath = os.path.expanduser("~/git/scqc")
    sys.path.append(gitpath)

    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)

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

    parser.add_argument('-q', '--query',
                        metavar='project_id',
                        type=str,
                        nargs='+',
                        default=None,
                        help='Perform standard query on supplied proj_id')

    parser.add_argument('-f', '--fasterq',
                        metavar='run_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')

    parser.add_argument('-p', '--download',
                        metavar='proj_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download (Runs) within project args with wget. e.g. SRR14584407')

    parser.add_argument('-m', '--metadata',
                        metavar='project_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download metadata for args. ')


    parser.add_argument('-i', '--impute',
                        metavar='proj_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='impute tech')


    parser.add_argument('-o', '--outfile',
                        metavar='outfile',
                        type=str,
                        default=None,
                        help='Outfile. ')

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

    if args.setup:
        s = setup(cp)
        # s.execute()

    if args.query is not None:
        q = Query(cp)
        for pid in args.query:
            q.execute(pid)

    if args.download is not None:
        # start a queue
        dq = Queue() 
        # loop through each project id 
        for proj_id in args.download:
            # get the runs associated with that project
            # run_ids = get_runs_for_project(cp, proj_id)
            
                # download the SRA binary file for the run
            fq = Download(cp)
            fq.execute(proj_id)
            # dq.put(fq.execute(proj_id))

        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('download', 'max_downloads'))

        # limit number of jobs 
        for n in range(md): 
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')

    if args.fasterq is not None:
        dq = Queue()
        for proj_id in args.fasterq:
            srr_ids = get_runs_for_project(cp, proj_id)
            for srr in srr_ids:
                fq = FasterqDump(cp, srr)
                dq.put(fq)
            logging.debug(f'created queue of {dq.qsize()} items')
            md = int(cp.get('analyze', 'max_jobs'))

        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')

    elif args.metadata is not None:
        for srp in args.metadata:
            df = query_project_metadata(srp)
            exps = list(df['Experiment'].unique())
            logging.debug(f"Got list of {len(exps)} experiments")
            for e in exps:
                print(e)

    if args.impute is not None:
        q = Impute(cp)
        for pid in args.impute:
            q.execute(pid)

