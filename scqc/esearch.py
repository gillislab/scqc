#!/usr/bin/env python
#
# Core utility to query esearch and create input project lists/ and 
# experiment lists that match search. 
#
#
#

import argparse
import logging
import io
import os

from configparser import ConfigParser

class SraSearch(object):
    
    def __init__(self):
        self.log = logging.getLogger()
        self.log.info('running search')


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

        subparsers = parser.add_subparsers(dest='subcommand',
                                           help='sub-command help.')

        parser_analysis = subparsers.add_parser('search',
                                                help='query daemon')

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
            s = SraSearch()
            s.run()
 

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


