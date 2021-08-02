#!/usr/bin/env python
#
# Core utility to query esearch and create input project lists/ and 
# experiment lists that match search. 
#

import argparse
import logging
import io
import json
import os
import requests
import time
import traceback

import xml.etree.ElementTree as et

from configparser import ConfigParser
from requests.exceptions import ChunkedEncodingError
from urllib import parse


class SraSearch(object):
    
    def __init__(self, config, species, strategy, textword):
        """
        config:    scqc configparser
        species:   List of Linnean species names. 
        strategy:  List of strategy strings. 
        textword:  List of text word strings. 
        
        sra_esearch=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra
        &datetype=pdat&mindate=2001 &maxdate=2020
        
        """
        self.log = logging.getLogger()
        self.log.info('running search')
        self.config = config
        self.esearch_rooturl=config.get('sra','sra_esearch')
        self.expid_file = os.path.expanduser(config.get('sra','expid_file'))
        self.batchsize = int(config.get('sra','uid_batchsize'))
        self.query_max = int(config.get('sra','query_max'))
        self.efetch_sleep = float(config.get('sra','query_sleep'))
        self.full_url = self.build_searchurl(species, strategy, textword)
        self.log.debug('SraSearch initted. ')


    def build_searchurl(self, species=["mus musculus"], 
                                 strategy=["rna seq"], 
                                 textword=["single cell","brain"]
                                 ):
        """
        Given 3 lists, build search string for URL with following logic:
        
        (sp1 OR sp2) AND (st1 or st2) AND (tw1 AND tw2) 

        """
        TAGS= { 'species': '[Organism]',
                'strategy' : '[Strategy]',
                'textword' : '[Text Word]'
              }
        tspec = []
        tstrat = []
        ttxt = []
        for s in species:
            tspec.append(f"{s}{TAGS['species']}")
        for r in strategy:
            tstrat.append(f"{r}{TAGS['strategy']}")
        for t in textword:
            ttxt.append(f"{t}{TAGS['textword']}")
        spc = ' OR '.join(tspec)
        strat = ' OR '.join(tstrat)
        txtwrd = ' AND '.join(ttxt)
        params = f"({spc}) AND ({strat}) AND ({txtwrd})"
        self.log.debug(f'uncoded params = {params}')
        encoded = parse.quote_plus(params)        
        furl = f"{self.esearch_rooturl}&term={encoded}"
        self.log.info(f'full esearch URL= {furl}')
        return furl    


    def run(self):
        self.log.info(f'running with url: {self.full_url}')
        uidlist = self.query_all_uids()
        self.log.debug(f'got {len(uidlist)} uids from search.')
        tuplist = []
        curid = 0
        while curid < len(uidlist):
            dolist = uidlist[curid:curid + self.batchsize]
            outtups = query_project_for_uidlist(self.config, dolist)
            for (expid, projid) in outtups:
                tuplist.append( (expid, projid)  )
                #if projid is not None and expid is not None :
                #    f.write(f'{expid} {projid}\n')
                #    f.flush()
                #else:
                #    self.log.warning('exp_id or proj_id is None. Ignoring... ')
            curid += self.batchsize
            self.log.debug(f'handled {self.batchsize} uids. Sleeping...')
            time.sleep(self.efetch_sleep)
        
        self.log.debug(f'got {len(tuplist)} tuples.')
        for (one, two ) in tuplist:
            print(f' {one} {two}')
        self.log.info(f'Done.')


    def query_all_uids(self):
        """
        Perform standard query with max return, loop over all entries. 
        retstart="+str(retstart)+"&retmax="+str(retmax)+"&retmode=json" 
    
        """
        #sra_esearch = config.get('sra', 'sra_esearch')
        #sra_efetch = config.get('sra', 'sra_efetch')
        #search_term = config.get('sra', 'search_term')
        query_max = 50000
        query_start = 0   
        log = logging.getLogger('sra')
        log.info('querying SRA...')
        alluids = []
        while True:
            try:
                url = f"{self.full_url}&retstart={query_start}&retmax={query_max}&retmode=json"
                log.debug(f"search url: {url}")
                r = requests.get(url)
                er = json.loads(r.content.decode('utf-8'))
                #log.debug(f"er: {er}")
                idlist = er['esearchresult']['idlist']
                log.debug(f"got idlist length={len(idlist)}")
                if len(idlist) > 0:
                    for id in idlist:
                        alluids.append(id)
                    query_start += query_max
                else:
                    break
            except Exception as ex:
                #log.error(f'problem with NCBI uid {id}')
                log.error(traceback.format_exc(None))
    
        log.debug(f'found {len(alluids)} uids.')
        return alluids


class SearchCLI(object):

    def runsearch(self):

        #FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
        #logging.basicConfig(format=FORMAT)

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

        parser.add_argument('-o', '--outfile',
                        metavar='outfile',
                        type=str,
                        default=None,
                        help='Outfile. ')
        
        subparsers = parser.add_subparsers(dest='subcommand',
                                           help='sub-command help.')

        parser_search = subparsers.add_parser('search',
                                                help='search command.')
    
        parser.add_argument('-u', '--uidquery',
                        metavar='uidfile',
                        type=str,
                        dest='uidquery',
                        required=False,
                        default=None,
                        help='Perform standard query on uids in file, print project_ids.')        

        parser_search.add_argument('-s','--species',
                            type=str,
                            dest='species',
                            default='mus musculus',
                            help='comma-separated list of lowercase Linnean names'
                            )

        parser_search.add_argument('-r','--strategy',
                            type=str,
                            dest='strategy',
                            default='rna seq',
                            help='comma-separated list of key words'
                            )

        parser_search.add_argument('-t','--textword',
                            type=str,
                            dest='textword',
                            default="'single cell','brain'",
                            help='comma-separated list of key words'
                            )
        
        args = parser.parse_args()

        # default to INFO
        logging.getLogger().setLevel(logging.INFO)

        if args.debug:
            logging.getLogger().setLevel(logging.DEBUG)
        if args.verbose:
            logging.getLogger().setLevel(logging.INFO)

        self.cp = ConfigParser()
        self.cp.read(os.path.expanduser(args.conffile))

        self.setuplogging(args.subcommand)           
        
        cs = self.get_configstr(self.cp)
        logging.debug(f"config: \n{cs} ")
        logging.debug(f"args: {args} ")

        if args.subcommand == 'search':
            sps = [x.strip() for x in args.species.split(',')]
            stg = [r.strip() for r in args.strategy.split(',')]
            txt = [t.strip() for t in args.textword.split(',')]
            
            s = SraSearch(self.cp, sps, stg, txt)
            s.run()

        if args.uidquery is not None:
            qfile = args.uidquery
            uidlist = readlist(os.path.expanduser(qfile))
            tuplist = []
            curid = 0
            batchsize = int(cp.get('sra','uid_batchsize'))
            efetch_sleep = float(cp.get('sra','query_sleep'))
            
            with open(args.outfile, 'w') as f:
                while curid < len(uidlist):
                    dolist = uidlist[curid:curid + batchsize]
                    outtups = query_project_for_uidlist(cp, dolist)
                    for (expid, projid) in outtups:
                        if projid is not None and expid is not None :
                            f.write(f'{expid} {projid}\n')
                            f.flush()
                        else:
                            logging.warning('exp_id or proj_id is None. Ignoring... ')
                    curid += batchsize
                time.sleep(efetch_sleep)

    def get_configstr(self, cp):
        with io.StringIO() as ss:
            cp.write(ss)
            ss.seek(0)  # rewind
            return ss.read()

    def setuplogging(self, command):
        """ 
        Setup logging 
 
        """
        self.log = logging.getLogger()

        #logfile = os.path.expanduser(self.cp.get(command, 'logfile'))
        #if logfile == 'syslog':
        #    logStream = logging.handlers.SysLogHandler('/dev/log')
        #elif logfile == 'stdout':
        logStream = logging.StreamHandler()
        #else:
        #    logStream = logging.FileHandler(filename=logfile, mode='a')    

        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        #formatter.converter = time.gmtime  # to convert timestamps to UTC
        logStream.setFormatter(formatter)
        self.log.addHandler(logStream)
        self.log.info('Logging initialized.')


    def run(self):
        self.runsearch() 




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


def query_project_for_uidlist(config, uidlist):
    """
    Non OOP version of the Query method
    """
    log = logging.getLogger('sra')
    sra_efetch = config.get('sra', 'sra_efetch')
    uids = ','.join(uidlist)
    url = f"{sra_efetch}&id={uids}"
    log.info(f"fetching url={url}")
    tries = 0
    while True:
        try:
            tuples = []
            r = requests.post(url)
            if r.status_code == 200:
                rd = r.content.decode()
                root = et.fromstring(rd)
                expkgs = root.findall('EXPERIMENT_PACKAGE')
                for exppkg in expkgs:
                    exp = exppkg.find('EXPERIMENT')
                    exp_id = exp.get('accession')
                    proj_id = exp.find('STUDY_REF').get('accession')
                    log.debug(f'exp_id: {exp_id} proj_id: {proj_id}')
                    tuples.append( (exp_id, proj_id) )
            time.sleep(0.5)
            return tuples
        
        except ChunkedEncodingError as cee:
            log.warning(f'got ChunkedEncodingError error for uidlist {uids}: {cee}')
            tries += 1
            if tries >= 3 and len(uidlist) == 1:
                log.warning(f'got ChunkedEncodingError for uid: {uidlist[0]} Giving up, returning None.')
                return None
            elif tries >= 3:
                log.warning(f'got too many ChunkedEncodingErrors for uidlist {uidlist}. Doing one-by-one...')
                tuples = query_project_for_uidlist_byone(config, uidlist)
                return tuples
            else:
                log.warning(f'got ChunkedEncodingError. Try {tries}')
            time.sleep(0.5)                
            
            
        except ConnectionError as ce:
            log.warning(f'got connection error for uidlist {uids}: {ce}')
            time.sleep(30)

        except Exception as e:
            log.warning(f'got another exception for uidlist {uids}: {e}  ')
            return None


def query_project_for_uidlist_byone(config, uidlist):
    """
    Fallback routine for troublesome encoding error with particular uids. 
    """
    log = logging.getLogger('sra')
    tuples = []
    for uid in uidlist:
        log.debug(f'handling single uid: {uid}')
        tuplist = query_project_for_uidlist(config, [uid])
        for tup in tuplist:
            if tup is not None:
                tuples.append(tup)
            else:
                log.warning(f'Problem querying during special handling. uid {uid}')
    log.warning(f'returning special tuplelist: {tuples}')
    return tuples
