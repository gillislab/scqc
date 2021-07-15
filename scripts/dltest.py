#!/usr/bin/env python
#
# Test script for testing download methods. 
#
#  https://sra-download.ncbi.nlm.nih.gov/traces/era22/ERR/ERR4822/ERR4822769
#  https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR4822769/ERR4822769


import argparse
import requests
import logging
import os
import sys
import time
import traceback

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import *
from scqc.sra import *






if __name__ == "__main__":
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

    parser.add_argument('-i', '--infile', 
                        metavar='infile', 
                        type=str, 
                        help='a flat file of run urls. ')   

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    urllist = readlist(args.infile)
    logging.info(f'got list of {len(urllist)} file_urls.')
    
