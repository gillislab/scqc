#!/usr/bin/env python
#

import argparse
import logging
import os 
import re
import sys
import tarfile

from configparser import ConfigParser


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *


PROJ_COLUMNS = ['proj_id', 'ext_ids', 'title', 'abstract', 'submission_id', 'data_source']

SAMP_COLUMNS = ['samp_id', 'ext_ids',  'taxon',
                'sciname', 'title', 'attributes', 'proj_id', 'submission_id', 'data_source']

EXP_COLUMNS = ['exp_id', 'ext_ids',  'strategy',
               'source', 'lcp', 'samp_id', 'proj_id', 'submission_id', 'data_source']

RUN_COLUMNS = ['run_id', 'ext_ids', 'tot_spots', 'tot_bases', 'run_size', 'publish_date',
               'taxon', 'organism', 'nreads',  'basecounts', 'file_url','file_size','exp_id', 'samp_id', 'proj_id', 
               'submission_id', 'data_source' ]

IMPUTE_COLUMNS = ['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon','batch', 'data_source']

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

class Impute(object):
    """
    Imputes sequencing technology for all runs under a project. 

    """
    def __init__(self, config):
        self.log = logging.getLogger('impute')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('impute', 'metadir'))
        self.resourcedir = os.path.expanduser(self.config.get('impute', 'resourcedir'))


    def execute(self, proj_id):
        """
        XXX if tech is not supported, will return an empty dataframe
        For proj_id:
            Impute technology where possible. 
            Put completed project ids into query-donelist.txt

        """
        self.log.info(f'handling proj_id {proj_id}')
        try:
            # read in experiment file
            expfile = f'{self.metadir}/experiments.tsv'
            edf = load_df(expfile)
            self.log.debug(f'opened experiments DF OK...')
            edf = edf[edf.proj_id == proj_id].reset_index(drop=True) # rename pdf -> edf 
            self.log.debug(f'got project-specific experiment df: \n{edf}')
            # auto impute technology  -  exp_id|tech

            idf = self.impute_tech_from_lcp(edf)              
            self.log.debug(f'got initial impute df: \n{idf}')
            # gather manually-curated smartseq tech
      
            # match run to tech
            runfile = f'{self.metadir}/runs.tsv'
            rdf = load_df(runfile)
            rdf = rdf[rdf.proj_id == proj_id].reset_index(drop=True)
            self.log.debug(f'got project-specific runs df: \n{rdf}')
            
            # impute 10x version
                
                        
            self.log.info(f'completed imputation for proj_id {proj_id}')
            #return (proj_id, proj_id, proj_id)          

        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            logging.error(traceback.format_exc(None))
            return (None, None, proj_id)
        
        
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


def setup(config):
    '''
    Builds directories in config file 
    Download the appropriate supplement data.
    Only needs to be done once.
    '''

    log = logging.getLogger('biccn')
    config = config
    # directories

    resourcedir = os.path.expanduser(config.get('setup', 'resourcedir'))
    outputdir = os.path.expanduser(config.get('setup', 'outputdir'))
    figuredir =os.path.expanduser(config.get('setup', 'figuredir'))

    dirs_to_make = [metadir, 
                    f'{cachedir}/nemo', 
                    tempdir, 
                    resourcedir, 
                    outputdir, 
                    figuredir]

    for direc in dirs_to_make:
        try:
            os.makedirs(direc)
        except FileExistsError:
            pass





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

    if args.impute is not None:
        q = Impute(cp)
        for pid in args.impute:
            q.execute(pid)
        


