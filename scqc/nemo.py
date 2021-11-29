#!/usr/bin/env python
#

import argparse
import logging
import os 
import re
import sys
import tarfile

from configparser import ConfigParser
from urllib.parse import urlparse

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *
from scqc.common import *

class Impute(object):
    """
    Imputes sequencing technology for all runs under a project. 

    """
    def __init__(self, config):
        self.log = logging.getLogger('impute')
        self.config = config
        self.cachedir = os.path.expanduser(self.config.get('impute','cachedir'))
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

            # get runs
            runfile = f'{self.metadir}/runs.tsv'
            rdf = load_df(runfile)
            rdf = rdf[rdf.proj_id == proj_id].reset_index(drop=True)
            self.log.debug(f'got project-specific runs df: \n{rdf}')
            
            # initial full def. default to whatever is in exp file (from LCP)
            idf = pd.merge(idf,rdf, left_on='exp_id',right_on='exp_id', how='left')
            idf['read1'] = ''
            idf['read2'] = ''
            idf['batch'] = '' 
            #idf['data_source'] = 'nemo'
                                 
            # impute 10x version
            idf = idf[IMPUTE_COLUMNS]
            #outdf.columns = IMPUTE_COLUMNS  # renames the columns from global 
            
            self.log.debug(f'initial full impute df: \n{idf}')            
            #outdf = self._known_tech(outdf)
            #idf.fillna(value='', inplace=True)
            
            idf = self.impute_tech_from_url(rdf, idf)
            self.log.debug(f'impute df after url inference: \n{idf}')            
            
            if len(idf) > 0:
                merge_write_df(idf, f'{self.metadir}/impute.tsv')  
            else :
                self.log.warn(f'Unable to predict tech for:{proj_id} ')
                return (None, None, proj_id)
            
            self.log.info(f'completed imputation for proj_id {proj_id}')
            return (proj_id, proj_id, proj_id)              

        except Exception as ex:
            self.log.error(f'problem with NCBI proj_id {proj_id}')
            logging.error(traceback.format_exc(None))
            return (None, None, proj_id)
        
        
    def impute_tech_from_lcp(self, edf):
        '''
        Take in experiment df for specific project.  
            Get unique library construction protocol (lcp) values. 
            For each value, loop over all keyword terms, identifying putative tech. 
            Fill in method column for all runs in temp DF.   
            Create impute DF
        
        '''
        self.log.debug(f'got edf: \n{edf}')
        edf.lcp = edf.lcp.fillna('').values

        ulcp = pd.DataFrame({"lcp": edf.lcp.unique()})
        ulcp['tech_version'] ='unknown'
        ulcp['techcount'] = 0
        # ulcp= ulcp.loc[1:,:].reset_index(drop=True)         
        # search for the keywords
        # doesn't play nice with NaN lcp - fill with a str

        for i in range(len(TECH_RES)):
            key = list(TECH_RES)[i]
            kw = TECH_RES[key]
            ulcp[key] = ulcp.lcp.str.contains(kw)
        
            ulcp.loc[ulcp[key],'tech_version'] = key
            ulcp.loc[ulcp[key],'techcount'] +=1
        
        ulcp.loc[ulcp['techcount'] > 1,'tech'] = 'multiple'
        dfout = ulcp[['lcp','tech_version']]
        self.log.debug(f'keyword hits: {ulcp.tech.values}')

        dfout = dfout.merge(edf[['exp_id','lcp']], on="lcp")
        return dfout[['exp_id','tech_version']]
    

    def impute_tech_from_url(self, rdf, idf):
        """
        need file_url and run_id from rdf. 
        set tech_version in idf 
        """
        #tdf = pd.merge(idf, rdf, left_on='run_id',right_on='run_id', how='left')
        idf.exp_id = idf.exp_id.astype('str')
        idf.samp_id = idf.samp_id.astype('str')
        tdf = pd.merge(idf, rdf, on=['run_id','proj_id','data_source','taxon','exp_id','samp_id'], how='left')
        
        for index, row in tdf.iterrows():
            runid = row['run_id']
            file_url = row['file_url']
            read1, read2, tech_version = self.process_tarfile(file_url) 
            tdf['read1'][index] = read1
            tdf['read2'][index] = read2
            tdf['tech_version'][index] = tech_version
        self.log.debug(f'got impute tech df:\n{tdf}')
        tdf = tdf[['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id','taxon' ,'batch', 'data_source']]
        return tdf    
                
                
    def process_tarfile(self, file_url):
        """
        E.g.
        file_url = 'https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOp/SM-DD44D_S32_E1-50.fastq.tar'
        ->    <cachedir>/nemo/<runid>.fastq.tar
        
        @return     read1,  read2,  tech_version -> 10xv1 | 10xv2 | 10xv3
        
        Strings, sequence length:

            SSv4     -> smargseq      >28
            10x_v1   -> 10xv1         24
            10x_v2   -> 10xv2         26
            10x_v3   -> 10xv3         28
        
        E.g.
        @D00201:CBDMVANXX170812:CBDMVANXX:7:1101:10000:57844 1:N:0:TAGCGCTCCTCTCTAT
        ATTAAAAGCAGTACTTAATTTGTGTTTCTCTGGCGCAAGTTTTTATCTTTG
        +
        BB/<BF<FBFFFFFFFFFFFBFFFFFFFFFFBFFFFFBFFFFFFFFFFFF#
        """
        #self.log.debug(f'handling file_url: {file_url}')
        o = urlparse(file_url)
        tbase = o.path.split('/')[-1:][0]
        tf = f"{self.cachedir}/nemo/{tbase}"
        #self.log.debug(f"handling tarfile {tf}")
        read1 = ''
        read2 = ''
        tech_version = self.scan_url_tech(o.path)
        to = tarfile.open(tf)
        subfiles = to.getnames()
        self.log.debug(f'got {len(subfiles)} files in tarfile...')
        if (tech_version == 'smartseq') and (len(subfiles) == 2):
            self.log.debug('handling smartseq...')
            read1 = subfiles[0]
            read2 = subfiles[1]
        elif ('10x' in tech_version) and (len(subfiles) < 4):
            self.log.debug('handling 10x...')
            for f in subfiles:
                #self.log.debug(f'handling subfile {f}')
                (err, out, rc) = peek_tarball(tf, f, 3)
                sline = out.split('\n')[1]
                self.log.info(f'got line 2:\n{sline}')
                if len(sline) > 23 and len(sline) < 29 :
                    # technical read -> read2
                    read2 = f
                elif len(sline) > 30:
                    # biological read -> read1 
                    read1 = f
        else:
            self.log.warning(f'wrong/no tech or unsupported number of files. ')                
        self.log.debug(f'tech_version = {tech_version} read1={read1} read2={read2}')
        return read1, read2, tech_version


    def scan_url_tech(self, path):
        """
          E.g. path  '/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOp/SM-DD44D_S32_E1-50.fastq.tar'
            SSv4     -> smartseq      >28
            10x_v1   -> 10xv1         24
            10x_v2   -> 10xv2         26
            10x_v3   -> 10xv3 
        """
        tech = ''
        for k in NEMO_URL_TECH_MAP.keys():
            if k in path:
                #self.log.debug(f"found {k} in {path}, returning tech={NEMO_URL_TECH_MAP[k]}")
                return NEMO_URL_TECH_MAP[k]
        return tech

def stage_in(config, cachedir, tempdir, runlist, force=True):
    """
    
    """
    runlen = len(runlist)
    logging.debug(f'handling runlist w/ {runlen} tarfiles...')

    for i, tf in enumerate(runlist):
        fpath = f"{cachedir}/nemo/{tf}.fastq.tar"
        to = tarfile.open(fpath)
        subfiles = to.getnames()
        for f in subfiles:
            if os.path.exists(f'{tempdir}/{f}') and not force:
                logging.debug(f'[{i}/{runlen}] path {tempdir}/{f} exists and force is not set. Skipping.')
            else:
                logging.debug(f'[{i}/{runlen}] extracting {f} to {tempdir}')
                to.extract(f, path=tempdir)
    logging.debug(f'done extracting files.')


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
        


