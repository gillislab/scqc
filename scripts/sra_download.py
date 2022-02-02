#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.sra import Query
from scqc.utils import *


#  Standard SRR download path:
#  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7117190/SRR7117190
# 

URLBASE='https://sra-pub-run-odp.s3.amazonaws.com/sra'

def download_run_sra(runid):

    srcurl=f'{URLBASE}/{runid}/{runid}'

    try:
        download_wget(srcurl, destpath=f'./{runid}.sra',rate = '50M')
    except KeyError:
        pass


# to do: include drivers for 10x and ss alignments
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

    parser.add_argument('-f', '--fasterq',
                        metavar='fastqrun',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')

    parser.add_argument('-r', '--runs',
                        metavar='run_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download run SRA file, e.g. SRR14584407')

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
    
    if args.runs is not None:
        for run_id in args.runs:
            logging.debug(f'Downloading SRA for run {run_id}')
            download_run_sra(run_id)
        logging.info('Done')
    

    
