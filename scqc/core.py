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
                self.todolist = self.readlist(self.todofile)
                self.donelist = self.readlist(self.donefile)
                if self.todolist is not None:
                    self.dolist = self.listdiff(self.todolist, self.donelist) 
                else:
                    self.dolist = []
                self.finished = self.execute(self.dolist)
                self.log.debug(f"got donelist of length {len(self.finished)}. writing...")
                self.writedone(self.finished)
                self.log.debug(f"done writing donelist: {self.donefile}. sleeping {self.sleep} ...")
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.log.error(traceback.format_exc(None))
            raise ex       

    def readlist(self, filepath):
        if filepath is not None:
            self.log.info(f'reading file: {filepath}')
            flist = []
            try:
                with open(filepath, 'r') as f:
                   flist = [line.strip() for line in  f]
                self.log.debug(f'got list with {len(flist)} items.')
                return flist 
            except:
                return []
        else:
            self.log.info('no file. return [].')
            return []

    def writedone(self, finishedlist):
        if self.donefile is not None:
            self.log.info('reading current done.')
            donelist = self.readlist(self.donefile)            
            self.log.info('adding finished.')
            alldone = self.listmerge(finishedlist, donelist)            
            self.log.info('writing donefile...')
            rootpath = os.path.dirname(self.donefile)
            basename = os.path.basename(self.donefile)
            try:
                (tfd, tfname) = tempfile.mkstemp(suffix=None, 
                                              prefix=f"{basename}.", 
                                              dir=f"{rootpath}/", 
                                              text=True)
                self.log.debug(f"made temp {tfname}")
                with os.fdopen(tfd, 'w') as f:
                    nlines = 0
                    for item in alldone:
                        f.write(f"{item}\n")
                        nlines += 1    
                os.rename(tfname, self.donefile)
                self.log.info(f"wrote {nlines} to {self.donefile}")
            except Exception as ex:
                self.log.error(traceback.format_exc(None))
                
            finally:
                pass

        else:
            self.log.info('no donefile defined.')


    def listdiff(self, list1, list2):
        self.log.debug(f"got list1: {list1} list2: {list2}")
        s1 = set(list1)
        s2 = set(list2)
        sd = s1 - s2
        dl = list(sd)
        dl.sort()
        self.log.debug(f"diff has length {len(dl)}")
        return dl

    def listmerge(self, list1, list2):
        self.log.debug(f"got list1: {list1} list2: {list2}")
        s1 = set(list1)
        s2 = set(list2)
        sd = s1 | s2
        dl = list(sd)
        dl.sort()
        self.log.debug(f"merged has length {len(dl)}")
        return dl
          
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
   

class Analysis(object):
    
    def __init__(self, config):
        self.log = logging.getLogger('analysis')
        self.log.info('analysis init...')
        self.config = config
        self.shutdown = False
        self.sleep = int(self.config.get('analysis','sleep')) 
        
    def run(self):
        self.log.info('analysis run...')
        try:
            while not self.shutdown:
                self.log.debug(f'analysis cycle. {self.sleep} seconds...')
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.stop()
            raise ex        
    
    def stop(self):
        self.log.info('stopping...')        


class Statistics(object):

    def __init__(self, config):
        self.log = logging.getLogger('statistics')
        self.log.info('statistics init...')
        self.config = config
        self.shutdown = False
        self.sleep = int(self.config.get('statistics','sleep')) 
        
    def run(self):
        self.log.info('statistics run...')
        try:
            while not self.shutdown:
                self.log.debug(f'statistics cycle. {self.sleep} seconds...')
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.stop()
            raise ex        
    
    def stop(self):
        self.log.info('stopping...')        



class CLI(object):
          
    def parseopts(self):   
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
    
    
    
            
            
            
