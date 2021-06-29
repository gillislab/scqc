
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
import pandas as pd
import numpy as np

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

PROJ_COLUMNS = ['proj_id', 'ext_ids', 'title', 'abstract', 'submission_id']

SAMP_COLUMNS = ['samp_id', 'ext_ids',  'taxon',
                'organism', 'title', 'attributes', 'proj_id', 'submission_id']

EXP_COLUMNS = ['exp_id', 'ext_ids',  'strategy',
               'source', 'lcp', 'samp_id', 'proj_id', 'submission_id']

RUN_COLUMNS = ['run_id', 'ext_ids', 'tot_spots', 'tot_bases', 'size', 'publish_date',
               'taxon', 'organism', 'nreads',  'basecounts', 'exp_id', 'samp_id', 'proj_id', 'submission_id']

IMPUTE_COLUMNS = ['proj_id','exp_id','samp_id','run_id', 'tech']
# run_id tech_version                read1                read2      exp_id   proj_id taxon
TECH_RES = {
    '10x'   : re.compile("10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'|Kit v1|PN-120233|10X V1", re.IGNORECASE),
    #'10xv1' : re.compile("", re.IGNORECASE),
    #'10xv2' : re.compile("v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2", re.IGNORECASE),
    #'10xv3' : re.compile("v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078", re.IGNORECASE),
    'smartseq'    : re.compile("Smart-Seq|SmartSeq|Picelli|SMART Seq", re.IGNORECASE),
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




def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()


class Impute(object):
    """
    Imputes sequencing technology for all runs under a project. 

    """
    def __init__(self, config):
        self.log = logging.getLogger('impute')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('impute', 'metadir'))
        # self.cachedir = os.path.expanduser(
        #     self.config.get('impute', 'cachedir'))
        # self.sra_esearch = self.config.get('sra', 'sra_esearch')
        # self.sra_efetch = self.config.get('sra', 'sra_efetch')
        # self.sleep = float(self.config.get('sra', 'sleep'))

    def execute(self, projectid):
        """
         For projectid:
            Impute technology where possible. 
            Put completed project ids into query-donelist.txt
            examples    - SRP114926 - contains both 10x and smartseq
                        - SRP122508 - contains just 10xv2. 192 runs, 10 exp, 10 samples
        """
        self.log.info(f'handling projectid {projectid}')
        try:
            # read in experiment file
            expfile = f'{self.metadir}/experiments.tsv'
            edf = pd.read_csv(expfile, sep='\t', index_col=0)
            self.log.debug(f'opened experiments DF OK...')
            edf = edf[edf.proj_id == projectid].reset_index(drop=True) # rename pdf -> edf 
            self.log.debug(f'got project-specific df: \n{edf}')
            # impute technology  -  exp_id|tech
            idf = self.impute_tech_from_lcp(edf)    
            self.log.debug(f'got initial imputed tech df: \n{idf}')

            # match run to tech
            runfile = f'{self.metadir}/runs.tsv'
            rdf = pd.read_csv(runfile, sep='\t', index_col=0)
            rdf = rdf[rdf.proj_id==projectid].reset_index(drop=True)
            # impute 10x version
            fdf = self.impute_10x_version(idf, rdf)
            self.log.debug(f'got imputed 10x version df: \n{fdf}')
            ssdat = self.parse_smartseq(idf, rdf)

            # save to disk
            merge_write_df(fdf, f'{self.metadir}/tech.tsv' )          
            self.log.info(f'completed imputation for proj_id {projectid}')
            return projectid          

        except Exception as ex:
            self.log.error(f'problem with NCBI projectid {projectid}')
            logging.error(traceback.format_exc(None))
            raise ex

    # old XXX
    def __impute_tech(self, df):
        ulcp = pd.DataFrame({"lcp": df.lcp.unique()})
    
        # search for the keywords
        
        for i in range(len(TECH_RES)):
            key = list(TECH_RES)[i]
            kw = TECH_RES[key]
            ulcp[key] = ulcp.lcp.str.contains(kw)


            ulcp[key] = ulcp.lcp.str.contains(
                kw, case=False, regex=True)
        self.log.debug(f'keyword hits: {ulcp.values}')
        
        tmpdf = ulcp.loc[:, list(TECH_RES.keys())]
        tmpdf = tmpdf.fillna(False)
        #unknownTechs = ulcp.loc[tmpdf.sum(axis=1) ==0,"LCP"]
        tmpdf["some10x"] = tmpdf.loc[:, "is10x"]
        # tmpdf["isMultiple"] = tmpdf.iloc[:,4:-1].sum(axis=1)  > 1
        tmpdf["smartseq"] = tmpdf.loc[:, "ss"] & ~(
            tmpdf.loc[:, "is10x"].astype('bool'))
        tmpdf["method"] = "Unknown"
        
        for tech in tmpdf.columns.values[5:-1]:
            tmpdf.loc[tmpdf.loc[:, tech], "method"] = tech
        
        ulcp['method'] = tmpdf.method
        ulcp = ulcp.loc[:, ["lcp", "method"]]
        df = df.merge(ulcp, on="lcp")
        return df

    # updated 6/29 jlee. 
    # TODO deal with multiple tech finds
    def impute_tech_from_lcp(self,df):
        '''
        Take in experiment df for specific project.  
            Get unique library construction protocol (lcp) values. 
            For each value, loop over all keyword terms, identifying putative tech. 
            Fill in method column for all runs in temp DF.   
            Create impute DF
        
        '''
        logging.debug(f'got df: \n{df}')
        df.lcp= df.lcp.fillna('None').values()

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



    # TODO include logging
    def impute_10x_version(self,idf,rdf):
        """
        For known 10x, get first part of fasta file and determine version.
        Only looks at the 10x portion if multiple techs used.
        Returns only 10x runs
        """

        # look at the first few lines of the fastq file.
        loglev = LOGLEVELS[self.log.getEffectiveLevel()]

        
        # require runs have rdf.nreads > 1. Otherwise, unable to impute
        rdf = rdf [rdf.nreads > 1]
        df = rdf.merge(idf, on = 'exp_id',how='left')

        # get all runs associated with the 10x inferred experiments
        runs = df.loc[ df.tech == '10x','run_id']
        if len(runs)  == 0:
            print('no runs imputable') 
            return None

        allRows = []
        # TODO multithread?
        for srrid in runs :

            cmd = ['fastq-dump',
                '--maxSpotId', '1',
                '--split-spot',
                '--stdout',
                '--log-level', f'{loglev}',
                srrid]      # don't assume that sra file exists. most likely wont

            # TODO suppress messages from cmd
            cp = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            dat = cp.stdout.read().decode("utf-8").split('\n')[:-1]

            # get the lengths of each read.
            lengths = {}
            it = 1
            for line in dat[0::4]:  # look at every 4 lines for the length of the read
                lengths[f'{srrid}_{it}.fastq'] = line.split("length=")[-1]
                it += 1

            l = [int(l) for l in lengths.values()]
            ind = l.index(max(l))

            # which is the longest? use as cDNA read
            read_bio = list(lengths.keys())[ind]
            
            # 10xv2 is typically 98 bp
            # 10xv3 is typically 91 bp
            tech = "10x"
            for i in range(len(lengths)):
                if l[i] == 24:
                    tech = "10xv1"
                    ind2 = i
                elif l[i] == 26:
                    ind2 = i
                    tech = "10xv2"
                elif l[i] == 28:
                    ind2 = i
                    tech = "10xv3"

            if tech == "10x":
                read_tech = None
            else:
                read_tech = list(lengths.keys())[ind2]
            allRows.append([srrid, tech, read_bio, read_tech])

        outdf = pd.DataFrame(allRows,columns=['run_id','tech_version','read1','read2'] )
        tax = rdf[['run_id','exp_id','proj_id','taxon']]

        outdf = outdf.merge(tax, on='run_id', how='inner') 
        
        # make sure to return in order!
        return outdf



    def parse_smartseq(self,idf,rdf):
        """
        For known smartseq, parse, build manifest
        """

        log = logging.getLogger('sra')
        # look at the first few lines of the fastq file.
        loglev = LOGLEVELS[log.getEffectiveLevel()]

        # get all runs associated with the smartseq inferred experiments
        df = rdf.merge(idf , on="exp_id", how = 'inner')    
        df = df[['run_id','taxon','nreads','exp_id','tech']]
        runs = df.loc[ df.tech == 'smartseq','run_id']
        
        outdf=pd.DataFrame({'run_id' : runs ,'tech_version':'smartseq'})

        # used to construct manifest and be consistent with 10x section
        outdf.loc[df.nreads > 1 ,'read1']  = outdf['run_id']+ '_1.fastq'
        outdf.loc[df.nreads > 1 ,'read2']  = outdf['run_id']+ '_2.fastq'
        outdf.loc[df.nreads == 1 ,'read1'] = outdf['run_id']+ '.fastq'
        outdf.loc[df.nreads == 1 ,'read1'] = '-'
        tax = rdf[['run_id','exp_id','proj_id','taxon']]

        outdf = outdf.merge(tax, on='run_id', how='inner') 
        
        # make sure to return in order!
        return outdf


    def impute_batch(sdf, rdf):
        '''
        This should only be run at the project level dataframes, but just in case,
        Splits by project id first and assigns a set batches to each 
        project id 
        '''
        
        cols = rdf.columns.tolist()
        new_rdf = pd.DataFrame(columns= cols.append('batch'))
        for proj , df in sdf.groupby('proj_id') :

            samp2batch = pd.DataFrame({ 'samp_id' : df.samp_id ,
                            'batch': pd.factorize(df['attributes'])[0]})
            tm_rdf = rdf.loc[ rdf.proj_id == proj ,:] 
            # merge these batches with the runs 
            tm_rdf=tm_rdf.merge(samp2batch, how ="left", on="samp_id")

            new_rdf = pd.concat([new_rdf, tm_rdf])

        return new_rdf
