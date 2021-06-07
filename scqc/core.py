#!/usr/bin/env python
#
#

import argparse
import fcntl
import io
import logging
import os
import tempfile
import time
import traceback

from configparser import ConfigParser
from queue import Queue

from scqc import sra
from scqc.utils import *

class Stage(object):
    '''
    Handles stage in pipeline. 
    Reads donelist. Reads todolist. Calculates diff. 
    Executes all for difflist. 
    Writes updated donelist. 
    
    '''
    
    def __init__(self, config, name):
        self.name = name
        self.log = logging.getLogger(self.name)
        self.log.info(f'{self.name} init...')
        self.config = config
        self.todofile = self.config.get(f'{self.name}','todofile')
        if self.todofile.lower().strip() == "none":
            self.todofile = None
        else:
            self.todofile = os.path.expanduser(self.todofile)         
        
        self.donefile = self.config.get(f'{self.name}','donefile')
        if self.donefile.lower().strip() == "none":
            self.donefile = None
        else:
            self.donefile = os.path.expanduser(self.donefile)    
        self.shutdown = False
        self.sleep = int(self.config.get(f'{self.name}','sleep'))
        self.outlist = [] 


    def run(self):
        self.log.info(f'{self.name} run...')
        try:
            while not self.shutdown:
                self.log.debug(f'{self.name} cycle will be {self.sleep} seconds...')
                self.todolist = readlist(self.todofile)
                self.donelist = readlist(self.donefile)
                if self.todolist is not None:
                    self.dolist = listdiff(self.todolist, self.donelist) 
                else:
                    self.dolist = []
                self.finished = self.execute(self.dolist)
                self.log.debug(f"got donelist of length {len(self.finished)}. writing...")
                
                if self.donefile is not None and len(self.finished) > 0:
                    logging.info('reading current done.')
                    donelist = readlist(self.donefile)            
                    logging.info('adding just finished.')
                    alldone = listmerge(self.finished, donelist)      
                    writelist( self.donefile, alldone)
                    self.log.debug(f"done writing donelist: {self.donefile}. sleeping {self.sleep} ...")
                else:
                    logging.info('donefile is None or no new processing. No output.')
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.log.error(traceback.format_exc(None))
            raise ex       


          
    def stop(self):
        self.log.info('stopping...')        

class Query(Stage):
    
    def __init__(self, config):
         super(Query, self).__init__(config, 'query')
         self.log.debug('super() ran. object initialized.') 

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        '''
        self.log.debug('executing...')
        self.log.info('ignoring dolist for query.')       
        sq = sra.Query(self.config)
        outlist = sq.execute()
        return outlist

# to do: pf.execute() followed by FasterqDump
class Download(Stage):

    def __init__(self, config):
        super(Download, self).__init__(config, 'download')
        self.log.debug('super() ran. object initialized.')
        self.max_downloads = int(self.config.get('download','max_downloads'))
        self.num_streams = int(self.config.get('download','num_streams')) 

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        '''
        self.log.debug(f'executing {self.name}')
        outlist = []
        dq = Queue()
        for runid in dolist:
            pf = sra.Prefetch(self.config, runid, outlist)
            dq.put(pf)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(self.config.get('sra','max_downloads'))
        for n in range(md):
            sra.Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')
        return outlist
   

class Analysis(Stage):

    def __init__(self, config):
        super(Download, self).__init__(config, 'analysis')
        self.log.debug('super() ran. object initialized.')
        
    def execute(self):
        pass

class Statistics(Stage):

    def __init__(self, config):
        super(Download, self).__init__(config, 'analysis')
        self.log.debug('super() ran. object initialized.')
        
    def execute(self):
        pass






class CLI(object):
          
    def parseopts(self):   
        
        
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
        logging.basicConfig(format=FORMAT)
        
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
    
        subparsers = parser.add_subparsers( dest='subcommand',
                                            help='sub-command help.')
    
        parser_analysis = subparsers.add_parser('query',
                            help='query daemon')
        
        parser_download = subparsers.add_parser('download',
                            help='download daemon')
    
        parser_analysis = subparsers.add_parser('analysis',
                            help='analysis daemon')

        parser_analysis = subparsers.add_parser('statistics',
                            help='statistics daemon')
    
        args= parser.parse_args()
        
        # default to INFO
        logging.getLogger().setLevel(logging.INFO)
        
        if args.debug:
            logging.getLogger().setLevel(logging.DEBUG)
        if args.verbose:
            logging.getLogger().setLevel(logging.INFO)
        
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile))
        cs = self.get_configstr(cp)
        logging.debug(f"config: {cs} ")
        
        logging.debug(f"args: {args} ")
    
        if args.subcommand == 'query':
            d = Query(cp) 
            d.run()
        
        if args.subcommand == 'download':
            d = Download(cp)
            d.run()

        if args.subcommand == 'analysis':
            d = Analysis(cp)
            d.run()    

        if args.subcommand == 'statistics':
            d = Statistics(cp)
            d.run()        
            
    def get_configstr(self, cp):
        with io.StringIO() as ss:
            cp.write(ss)
            ss.seek(0) # rewind
            return ss.read()
    
    def run(self):
        self.parseopts()
    
    
    
            
            
            
