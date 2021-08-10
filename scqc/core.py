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

from scqc import sra, star
from scqc.utils import *


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp


class Stage(object):
    '''
    Handles stage in pipeline. 
    Reads donelist. Reads todolist. Calculates diff. 
    Executes all for difflist. 
    Writes updated donelist, seenlist. 

    .run() should not be overriden. 
    .execute() can be overriden, but doesn't need to be for normal case. 


    Multi-server setup. Hash on proj_id list. Need stable hash.
    https://death.andgravity.com/stable-hashing
    

    '''

    def __init__(self, config, name):
        self.name = name
        self.log = logging.getLogger(self.name)
        self.log.info(f'{self.name} init...')
        self.config = config
        self.todofile = self.config.get(f'{self.name}', 'todofile')
        self.seenfile = self.config.get(f'{self.name}', 'seenfile')
        if self.todofile.lower().strip() == "none":
            self.todofile = None
        else:
            self.todofile = os.path.expanduser(self.todofile)

        self.donefile = self.config.get(f'{self.name}', 'donefile')
        if self.donefile.lower().strip() == "none":
            self.donefile = None
        else:
            self.donefile = os.path.expanduser(self.donefile)
        self.seenfile = self.config.get(f'{self.name}', 'seenfile')
        if self.seenfile.lower().strip() == "none":
            self.seenfile = None
        else:
            self.seenfile = os.path.expanduser(self.seenfile)
        self.shutdown = False
        self.sleep = int(self.config.get(f'{self.name}', 'sleep'))
        self.batchsize = int(self.config.get(f'{self.name}', 'batchsize'))
        self.batchsleep = float(self.config.get(f'{self.name}', 'batchsleep'))
        self.ncycles = int(self.config.get(f'{self.name}', 'ncycles'))
        self.num_servers = int(self.config.get(f'{self.name}','num_servers'))
        self.server_index = int(self.config.get(f'{self.name}','server_idx'))
        self.outlist = []

    def run(self):
        self.log.info(f'{self.name} run...')
        cycles = 0
        try:
            while not self.shutdown:
                self.log.debug(
                    f'{self.name} cycle will be {self.sleep} seconds...')
                self.todolist = readlist(self.todofile)
                self.donelist = readlist(self.donefile)
                if self.todolist is not None:
                    self.dolist = listdiff(self.todolist, self.donelist)
                else:
                    self.dolist = []
                # if multi-server setup, filter todolist...
                if self.num_servers > 1:
                    self.dolist = modulo_filter(self.dolist, self.num_servers , self.server_index)
                    
                # cut into batches and do each separately, updating donelist. 
                logging.debug(f'dolist len={len(self.dolist)}')
                curid = 0
                while curid < len(self.dolist):
                    dobatch = self.dolist[curid:curid + self.batchsize]
                    logging.debug(f'made dobatch length={len(dobatch)}')
                    logging.debug(f'made dobatch: {dobatch}')
                    (finished, seen) = self.execute(dobatch)
                    try:
                        finished.remove(None)
                    except:
                        self.log.warn('Got None in finished list from an execute. Removed.')
                    self.log.debug(f"got finished list len={len(finished)}. writing...")
                    
                    if self.donefile is not None:
                        logging.info('reading current done.')
                        donelist = readlist(self.donefile)
                        logging.info('adding just finished.')
                        alldone = listmerge(finished, donelist)
                        writelist(self.donefile, alldone)
                        self.log.debug(
                            f"done writing donelist: {self.donefile}.")
                    else:
                        logging.info(
                            'donefile is None or no new processing. No output.')
                    
                    if self.seenfile is not None and len(seen) > 0:
                        logging.info('reading current seen.')
                        seenlist = readlist(self.seenfile)
                        logging.info('adding just seen.')
                        allseen = listmerge(seen, seenlist)
                        writelist(self.seenfile, allseen)
                        self.log.debug(
                            f"done writing seenlist: {self.seenfile}. sleeping {self.batchsleep} ...")
                    else:
                        logging.info(
                            'seenfile is None or no new processing. No output.')
                        
                    curid += self.batchsize
                    time.sleep(self.batchsleep)
                cycles += 1
                if cycles >= self.ncycles:
                    self.shutdown = True

                # overall stage sleep
                if not self.shutdown:
                    logging.info(f'done with all batches. Sleeping for stage. {self.sleep} sec...')
                    time.sleep(self.sleep)

        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')

        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.log.error(traceback.format_exc(None))
            raise ex
        logging.info(f'Shutdown set. Exiting {self.name}')


    def stop(self):
        self.log.info('stopping...')


class Query(Stage):
    """
    Stage takes in list of NCBI project ids. 
    Collects metadata on projects, samples, experiments, and runs. Stores in DFs. 
    Outputs complete project ids. 
    """

    def __init__(self, config):
        super(Query, self).__init__(config, 'query')
        self.log.debug('super() ran. object initialized.')

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        '''
        self.log.debug(f'performing custom execute for {self.name}')
        self.log.debug(f'got dolist len={len(dolist)}. executing...')
        outlist = []
        seenlist = []
        for projectid in dolist:
            self.log.debug(f'handling id {projectid}...')
            try:
                sq = sra.Query(self.config)
                (out, seen) = sq.execute(projectid)
                self.log.debug(f'done with {projectid}')
                if out is not None:
                    outlist.append(out)
                if seenlist is not None:
                    seenlist.append(seen)
            except Exception as ex:
                self.log.warning(f"exception raised during project query: {projectid}")
                self.log.error(traceback.format_exc(None))
        self.log.debug(f"returning outlist len={len(outlist)} seenlist len={len(seenlist)}")
        return (outlist, seenlist)


    def setup(self):
        sra.setup(self.config)

class Impute(Stage):
    """
    Stage takes in list of NCBI project ids. 
    Examines Library Construction Protocol, and where needed downloads first X kilobytes of run files to guess
    library technology. 
    
    Outputs complete project ids. 
    
    """
    def __init__(self, config):
        super(Impute, self).__init__(config, 'impute')
        self.log.debug('super() ran. object initialized.')

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        '''
        self.log.debug(f'performing custom execute for {self.name}')
        self.log.debug(f'got dolist len={len(dolist)}. executing...')
        outlist = []
        seenlist = []
        for projectid in dolist:
            self.log.debug(f'handling id {projectid}...')
            try:
                si = sra.Impute(self.config)
                (out, seen) = si.execute(projectid)
                self.log.debug(f'done with {projectid}')
                if out is not None:
                    outlist.append(out)
                if seenlist is not None:
                    seenlist.append(seen)
            except Exception as ex:
                self.log.warning(f"exception raised during project query: {projectid}")
                self.log.error(traceback.format_exc(None))
        self.log.debug(f"returning outlist len={len(outlist)} seenlist len={len(seenlist)}")
        return (outlist, seenlist)


    def setup(self):
        sra.setup(self.config)


class Download(Stage):

    def __init__(self, config):
        super(Download, self).__init__(config, 'download')
        self.log.debug('super() ran. object initialized.')
        self.max_downloads = int(self.config.get('download', 'max_downloads'))
        self.num_streams = int(self.config.get('download', 'num_streams'))

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        
        '''
        self.log.debug(f'performing custom execute for {self.name}')
        self.log.debug(f'got dolist len={len(dolist)}. executing...')
        outlist = []
        seenlist = []
        for projectid in dolist:
            self.log.debug(f'handling id {projectid}...')
            try:
                sd = sra.Download(self.config)
                (done, seen) = sd.execute(projectid)
                self.log.debug(f'done with {projectid}')
                if done is not None:
                    outlist.append(done)
                if seenlist is not None:
                    seenlist.append(seen)
            except Exception as ex:
                self.log.warning(f"exception raised during project query: {projectid}")
                self.log.error(traceback.format_exc(None))
        self.log.debug(f"returning outlist len={len(outlist)} seenlist len={len(seenlist)}")
        return (outlist, seenlist)


    def setup(self):
        sra.setup(self.config)


class Analyze(Stage):

    def __init__(self, config):
        super(Analyze, self).__init__(config, 'analyze')
        self.log.debug('super() ran. object initialized.')

    def execute(self, dolist):
        '''
        Perform one run for stage.  
        '''
        self.log.debug(f'performing custom execute for {self.name}')
        self.log.debug(f'got dolist len={len(dolist)}. executing...')
        outlist = []
        seenlist = []
        for projectid in dolist:
            self.log.debug(f'handling id {projectid}...')
            try:
                ar = star.AlignReads(self.config)
                (out, seen) = ar.execute(projectid)
                self.log.debug(f'done with {projectid}')
                if out is not None:
                    outlist.append(out)
                if seenlist is not None:
                    seenlist.append(seen)
            except Exception as ex:
                self.log.warning(f"exception raised during project query: {projectid}")
                self.log.error(traceback.format_exc(None))
        self.log.debug(f"returning outlist len={len(outlist)} seenlist len={len(seenlist)}")
        return (outlist, seenlist)


    def setup(self):
        star.setup(self.config)


class Statistics(Stage):

    def __init__(self, config):
        super(Statistics, self).__init__(config, 'analysis')
        self.log.debug('super() ran. object initialized.')

    def execute(self):
        self.log.debug(f'performing custom execute for {self.name}')
        self.log.debug(f'got dolist len={len(dolist)}. executing...')
        outlist = []
        seenlist = []
        for projectid in dolist:
            self.log.debug(f'handling id {projectid}...')
            try:
                st = statistics.Statistics(self.config)
                (out, seen) = st.execute(projectid)
                self.log.debug(f'done with {projectid}')
                if out is not None:
                    outlist.append(out)
                if seenlist is not None:
                    seenlist.append(seen)
            except Exception as ex:
                self.log.warning(f"exception raised during project query: {projectid}")
                self.log.error(traceback.format_exc(None))
        self.log.debug(f"returning outlist len={len(outlist)} seenlist len={len(seenlist)}")
        return (outlist, seenlist)


    def setup(self):
        pass




class CLI(object):

    def parseopts(self):

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

        parser.add_argument('-N', '--nocleanup',
                            action="store_true",
                            dest='nocleanup',
                            help='do not remove temp files [False]')

        parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default='~/git/scqc/etc/scqc.conf',
                            help='Config file path [~/git/scqc/etc/scqc.conf]')

        parser.add_argument('-s', '--setup',
                            action="store_true",
                            dest='setup',
                            help='perform setup for chosen daemon and exit...'
                            )
        
        parser.add_argument('-n','--ncycles',
                            action='store',
                            dest='ncycles',
                            default=None,
                            help='halt after N cycles'
                            )


        subparsers = parser.add_subparsers(dest='subcommand',
                                           help='sub-command help.')

        parser_analysis = subparsers.add_parser('query',
                                                help='query daemon')

        parser_analysis = subparsers.add_parser('impute',
                                                help='impute daemon')

        parser_download = subparsers.add_parser('download',
                                                help='download daemon')

        parser_analysis = subparsers.add_parser('analyze',
                                                help='analysis daemon')

        parser_analysis = subparsers.add_parser('statistics',
                                                help='statistics daemon')

        args = parser.parse_args()

        # default to INFO
        logging.getLogger().setLevel(logging.INFO)

        if args.debug:
            logging.getLogger().setLevel(logging.DEBUG)
        if args.verbose:
            logging.getLogger().setLevel(logging.INFO)

        self.cp = ConfigParser()
        self.cp.read(os.path.expanduser(args.conffile))
        
        if args.ncycles is not None:
            self.cp.set('DEFAULT','ncycles',str(int(args.ncycles)))
            self.cp.set(args.subcommand, 'logfile', 'stdout')

        if args.nocleanup:
            self.cp.set('DEFAULT','nocleanup', "True")
        
        if args.setup:
            self.cp.set(args.subcommand, 'logfile', 'stdout')
        
        self.setuplogging(args.subcommand)           
        
        cs = self.get_configstr(self.cp)
        logging.debug(f"config: \n{cs} ")
        logging.debug(f"args: {args} ")

        if args.subcommand == 'query':
            d = Query(self.cp)
            if args.setup:
                d.setup()
            else:
                d.run()

        if args.subcommand == 'impute':
            d = Impute(self.cp)
            if args.setup:
                d.setup()
            else:
                d.run()

        if args.subcommand == 'download':
            d = Download(self.cp)
            if args.setup:
                d.setup()
            else:
                d.run()

        if args.subcommand == 'analyze':
            d = Analyze(self.cp)
            if args.setup:
                d.setup()
            else:
                d.run()

        if args.subcommand == 'statistics':
            d = Statistics(self.cp)
            if args.setup:
                d.setup()
            else:
                d.run()


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

        logfile = os.path.expanduser(self.cp.get(command, 'logfile'))
        if logfile == 'syslog':
            logStream = logging.handlers.SysLogHandler('/dev/log')
        elif logfile == 'stdout':
            logStream = logging.StreamHandler()
        else:
            logStream = logging.FileHandler(filename=logfile, mode='a')    

        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        #formatter.converter = time.gmtime  # to convert timestamps to UTC
        logStream.setFormatter(formatter)
        self.log.addHandler(logStream)
        self.log.info('Logging initialized.')


    def run(self):
        self.parseopts()
