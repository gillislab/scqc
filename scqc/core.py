#!/usr/bin/env python
#
#

import argparse
import io
import logging
import os
import time

from configparser import ConfigParser

class Query(object):
    '''
    inlist=None
    max_downloads=5
    num_streams=1
    donelist=~/git/scqc/test/download-donelist.txt
    
    '''
    def __init__(self, config):
        self.log = logging.getLogger('query')
        self.log.info('query init...')
        self.config = config
        self.inlist = self.config.get('query','inlist')
        if self.inlist == "None":
            self.inlist = None
        self.shutdown = False
        self.sleep = int(self.config.get('query','sleep')) 
        
    def run(self):
        self.log.info('query run...')
        try:
            while not self.shutdown:
                self.log.debug(f'query cycle. {self.sleep} seconds...')
                self.execute()
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.stop()
            raise ex        
    
    def execute(self):
        '''
        Perform one run. 
        
        '''
        self.log.debug('executing...')
        
    
    
    
    def stop(self):
        self.log.info('stopping...')        


        
class Download(object):
    '''
    inlist=~/git/scqc/test/query-donelist.txt
    max_downloads=5
    num_streams=1
    donelist=~/git/scqc/test/download-donelist.txt
    
    '''
    def __init__(self, config):
        self.log = logging.getLogger('download')
        self.log.info('downloader init...')
        self.config = config
        self.shutdown = False
        self.sleep = int(self.config.get('download','sleep')) 
        self.idlist = os.path.expanduser(self.config.get('download','inlist'))
        self.max_downloads = int(self.config.get('download','max_downloads'))
        self.num_streams = int(self.config.get('download','num_streams'))
        self.donelist = os.path.expanduser(self.config.get('download','donelist'))


        
    def run(self):
        self.log.info('downloader run...')
        try:
            while not self.shutdown:
                self.log.debug(f'download cycle. {self.sleep} seconds...')
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.stop()
            raise ex        
    
    def stop(self):
        self.log.info('stopping...')        



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
    
    
    
            
            
            
