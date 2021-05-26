#!/usr/bin/env python
#
# http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=
#
# Could use  SRR14584407 SRR14584408 in example..

import argparse
import io
import itertools
import json
import logging
import os
from queue import Queue

import requests
import subprocess
import sys
import time 

from configparser import ConfigParser
from threading import Thread
from queue import Queue, Empty

import xml.etree.ElementTree as et
import pandas as pd

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


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp

def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0) # rewind
        return ss.read()

class Worker(Thread):
    def __init__(self, q):
        self.q = q
        super(Worker, self).__init__()

    def run(self):
        # Race condition, just try!
        while True:
            try:
                job = self.q.get_nowait()
                job.execute()
                self.q.task_done()
            except Empty:
                return

class Query(object):
    
    def __init__(self, config):
        self.log = logging.getLogger('sra')
        self.config = config
        self.metadir = os.path.expanduser(self.config.get('query','metadir'))
        self.cachedir = os.path.expanduser(self.config.get('query','cachedir'))
        self.sra_esearch = self.config.get('sra','sra_esearch')
        self.sra_efetch = self.config.get('sra','sra_efetch')
        self.search_term=self.config.get('sra','search_term')
        self.query_max=self.config.get('sra','query_max')
        self.uidfile=os.path.expanduser(self.config.get('sra','uidfile'))


    def execute(self):
        self.log.info('querying SRA...')
        url =f"{self.sra_esearch}&term={self.search_term}&retmax={self.query_max}&retmode=json" 
        self.log.debug(f"search url: {url}")
        r = requests.get(url)
        er = json.loads(r.content.decode('utf-8'))
        logging.debug(f"er: {er}")
        idlist = er['esearchresult']['idlist']
        logging.debug(f"got idlist: {idlist}")
        # filter ids by already done. 
        donelist = readlist(self.uidfile)
        idlist = listdiff(idlist, donelist)        
      
        if len(idlist) > 0:
            allrows = []
            doneids = []
            for id in idlist:
                try:
                    url=f"{self.sra_efetch}&id={id}"
                    self.log.debug(f"fetch url={url}")
                    r = requests.post(url)
                    rd = r.content.decode()
                    #logging.debug(f"data for id={id}: {rd}")
                    rows = self._parse_experiment_pkg(rd)
                    allrows = itertools.chain(allrows, rows)
                    doneids.append(id)
                except Exception as ex:
                    self.log.error(f'problem with NCBI uid {id}')
                    logging.error(traceback.format_exc(None))
                    time.sleep(1)
            
            newdone = listmerge(donelist, doneids)
            self.log.info(f'updating uid done list...')
            writelist(self.uidfile, newdone)
            
            df = pd.DataFrame(allrows , columns = ["project","experiment","submission", "runs","date","taxon_id", "organism","lcp","title","abstract" ])
            df = df.fillna(value = "")  # fill None with empty strings. 
            df["status"] = "uidfetched"
            
            df = self._impute_tech(df)
            
            filepath = f"{self.metadir}/all_metadata.tsv"
            logging.info(f"saving metadata df to {filepath}")
            df.to_csv( filepath,
                      sep="\t", 
                      mode = 'a', 
                      index=False, 
                      header=not os.path.exists(filepath))
    
            sl = list(df.runs)
            sralist = list(itertools.chain.from_iterable(sl))
            return sralist
        else:
            self.log.info('no new uids to process...')
            return []
        
    
    def _parse_experiment_pkg(self, xmlstr):
        root = et.fromstring(xmlstr)
        self.log.debug(f"root={root}")
        rows = []
        for exp in root.iter("EXPERIMENT_PACKAGE"):
            for lcp in exp[0].iter("LIBRARY_CONSTRUCTION_PROTOCOL") :
                lcp = lcp.text
            SRXs = exp[0].get('accession')
            SRAs = exp[1].get('accession')
            SRPs = exp[3].get('accession')
            # title = exp[3][1][0].text
            abstract =exp[3][1][2].text
            SRRs = []
            date=[]
            taxon=[]
            orgsm =[]
            for study in exp[3].iter("STUDY_TITLE") :
                title = study.text 
            for study in exp[3].iter("STUDY_ABSTRACT") :
                abstract = study.text 
    
            for run in exp.iter("RUN") : 
                SRRs.append(run.attrib['accession'])
                date.append(run.attrib['published'])
                for mem in run.iter("Member") : 
                    taxon.append(mem.attrib['tax_id'])
                    orgsm.append(mem.attrib['organism'])       
            row = [SRPs, SRXs, SRAs, SRRs, date, taxon, orgsm, lcp, title,abstract]
            self.log.debug(f'got SRRs: {SRRs}')
            rows.append(row)
        return rows

    def _impute_tech(self, df):
        '''
        Take in df. 
            Get unique library construciton protocol (lcp) values. 
            For each value, loop over all keyword terms, identifying tech. 
            Fill in method  column for all runs in temp DF. Merge DF.
        Return filled in DF.  
        
        '''
        keywords = {
            "is10x"    : "10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'",
            "v3"       : "v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078",
            "v2"       : "v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2",
            "v1"       : "Kit v1|PN-120233|10X V1",
            "ss"       : "Smart-Seq|SmartSeq|Picelli|SMART Seq",
            "smarter"  : "SMARTer",
            "dropseq"  : "Cell 161, 1202-1214|Macosko|dropseq|drop-seq",
            "celseq"   : "CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq",
            "sortseq"  : "Sort-seq|Sortseq|Sort seq",
            "seqwell"  : "Seq-Well|seqwell",
            "biorad"   : "Bio-Rad|ddSeq",
            "indrops"  : "inDrop|Klein|Zilionis",
            "marsseq2" : "MARS-seq|MARSseq|Jaitin et al|jaitin",
            "tang"     : "Tang",
            # "TruSeq":"TruSeq",
            "splitseq" : "SPLiT-seq",
            "microwellseq":"Microwell-seq"
        }
    
        ulcp = pd.DataFrame({"lcp" : df.lcp.unique() })
        
        # search for the keywords
        for i in range(len(keywords)) :
            key= list(keywords)[i]
            kw = keywords[key] 
            ulcp[key] = ulcp.lcp.str.lower().str.contains(kw, case=False, regex=True)
        self.log.debug(f'keyword hits: {ulcp.values}' )
        
        tmpdf = ulcp.loc[:,list(keywords.keys())]
        tmpdf = tmpdf.fillna(False)   
        #unknownTechs = ulcp.loc[tmpdf.sum(axis=1) ==0,"LCP"]
        tmpdf["issome10x"]  = tmpdf.loc[:,"is10x"]
        # tmpdf["isMultiple"] = tmpdf.iloc[:,4:-1].sum(axis=1)  > 1 
        tmpdf["isss"]    = tmpdf.loc[:,"ss"] & ~(tmpdf.loc[:,"is10x"].astype('bool'))
        tmpdf["method"] = "Unknown"

        for tech in tmpdf.columns.values[5:-1] :
            tmpdf.loc[tmpdf.loc[:,tech],"method"]  = tech

        ulcp['method'] = tmpdf.method
        ulcp = ulcp.loc[:,["lcp","method"]]
        df = df.merge(ulcp , on = "lcp")    
        return df    

        


class Prefetch(object):
    '''
        Simple wrapper for NCBI prefetch
    Usage: prefetch [ options ] [ accessions(s)... ]
    Parameters:  
        accessions(s)    list of accessions to process
    Options:
      -T|--type <file-type>            Specify file type to download. Default: sra
      -N|--min-size <size>             Minimum file size to download in KB
                                        (inclusive).
      -X|--max-size <size>             Maximum file size to download in KB
                                         (exclusive). Default: 20G
      -f|--force <no|yes|all|ALL>      Force object download - one of: no, yes,
                                         all, ALL. no [default]: skip download if
                                         the object if found and complete; yes:
                                         download it even if it is found and is
                                         complete; all: ignore lock files (stale
                                         locks or it is being downloaded by
                                         another process - use at your own
                                         risk!); ALL: ignore lock files, restart
                                         download from beginning
      -p|--progress                    Show progress
      -r|--resume <yes|no>             Resume partial downloads - one of: no, yes
                                         [default]
      -C|--verify <yes|no>             Verify after download - one of: no, yes
                                         [default]
      -c|--check-all                   Double-check all refseqs
      -o|--output-file <file>          Write file to <file> when downloading
                                         single file
      -O|--output-directory <directory>
                                       Save files to <directory>/
         --ngc <path>                  <path> to ngc file
         --perm <path>                 <path> to permission file
         --location <location>         location in cloud
         --cart <path>                 <path> to cart file
      -V|--version                     Display the version of the program
      -v|--verbose                     Increase the verbosity of the program
                                         status messages. Use multiple times for
                                         more verbosity.
      -L|--log-level <level>           Logging level as number or enum string.
                                         One of
                                         (fatal|sys|int|err|warn|info|debug) or
                                         (0-6) Current/default is warn
         --option-file file            Read more options and parameters from the
                                         file.
      -h|--help                        print this message


    '''
    def __init__(self, config, runid, outlist):
        self.log = logging.getLogger('sra')
        self.config = config
        self.runid = runid
        self.outlist = outlist
        self.sracache = os.path.expanduser(self.config.get('sra','cachedir'))
        self.log.debug(f'prefetch id {runid}')

    def execute(self):
        self.log.debug(f'prefetch id {self.runid}')
        loglev = LOGLEVELS[ self.log.getEffectiveLevel()]
        cmd = [     'prefetch', 
                    '-O' , f'{self.sracache}/', 
                    '--log-level', f'{loglev}', 
                    f'{self.runid}' ]
        cmdstr = " ".join(cmd)
        logging.debug(f"prefetch command: {cmd} running...")        
        cp = subprocess.run(cmd)
        logging.debug(f"Ran cmd='{cmd}' returncode={cp.returncode} {type(cp.returncode)} " )
        if str(cp.returncode) == "0":
            self.outlist.append(self.runid)
        

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
        -S|--split-files                 write reads into different files
        -v|--verbose                     Increase the verbosity of the program
                                         status messages. Use multiple times for
                                         more verbosity.
        
    '''

    def __init__(self, config, srrid):
        self.log = logging.getLogger('sra')
        self.id = srrid
        self.log.debug(f'downloading id {srrid}')


    def execute(self):
        self.log.debug(f'downloading id {self.id}')

 
def get_run_metadata(sraproject):
    '''
    
    '''
    url="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi"

    headers = {
        "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding":"gzip,deflate,sdch",
        "Accept-Language":"en-US,en;q=0.8",
        "Cache-Control":"no-cache",
        "Connection":"keep-alive",
        "DNT":"1",
        "Host":"trace.ncbi.nlm.nih.gov",
        "Origin":"http://trace.ncbi.nlm.nih.gov",
        "Pragma":"no-cache",
        "Referer":"http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?db=sra",
        "User-Agent":"Mozilla/5.0 (iPhone; CPU iPhone OS 6_0 like Mac OS X) AppleWebKit/536.26 (KHTML, like Gecko) Version/6.0 Mobile/10A5376e Safari/8536.25"}
    
    payload = {
        "db":"sra",
        "rettype":"runinfo",
        "save":"efetch",
        "term": sraproject }
    
    r = requests.put(url, data=payload, headers=headers, stream=True) 
    with io.BytesIO(r.content) as imf:
        df = pandas.read_csv(imf)
    return df


if __name__ == "__main__":

    gitpath=os.path.expanduser("~/git/scqc")
    sys.path.append(gitpath)


    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
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
      
    parser.add_argument('-q','--query',
                        action="store_true", 
                        dest='query', 
                        help='Perform standard query')    
    
    parser.add_argument('-f','--fasterq', 
                        metavar='fasterq', 
                        type=str,
                        nargs='+',
                        required=False,
                        default=None, 
                        help='Download args with fasterq-dump. e.g. SRR14584407')  

    parser.add_argument('-p','--prefetch', 
                        metavar='prefetch', 
                        type=str,
                        nargs='+',
                        required=False,
                        default=None, 
                        help='Download args with prefetch. e.g. SRR14584407') 

    
    parser.add_argument('-m','--metadata',
                        metavar='metadata', 
                        type=str,
                        nargs='+',
                        required=False,
                        default=None, 
                        help='Download metadata for args. ')  
    
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
    
    cp = get_default_config()
    cs = get_configstr(cp)
    
    logging.debug(f"got config: {cs}")

    if args.query:
        q = Query(cp)
        q.execute()

    if args.fasterq is not None:
        dq = Queue()
        for srr in args.fasterq:
            fq = FasterqDump(cp, srr)
            dq.put(fq)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('sra','max_downloads'))
        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')

    if args.prefetch is not None:
        dq = Queue()
        for srr in args.prefetch:
            fq = Prefetch(cp, srr)
            dq.put(fq)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(cp.get('sra','max_downloads'))
        for n in range(md):
            Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')


    elif args.metadata is not None:
        for srr in args.metadata:
            df = get_run_metadata(srr)
            logging.debug(f"Got list of {len(df)} runs")
            runs = list(df['Run'])
            logging.info(f"Runlist: {runs}")
