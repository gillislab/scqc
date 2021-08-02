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


def download_wget(srcurl, destpath, finalname=None, overwrite=True, decompress=True, rate='1M'):
    """
    GNU Wget 1.20.1, a non-interactive network retriever.
    Usage: wget [OPTION]... [URL]...
    
    Startup:
      -V,  --version                   display the version of Wget and exit
      -h,  --help                      print this help
      -v,  --verbose                   be verbose (this is the default)
      -nv, --no-verbose                turn off verboseness, without being quiet
           --report-speed=TYPE         output bandwidth as TYPE.  TYPE can be bits
      -t,  --tries=NUMBER              set number of retries to NUMBER (0 unlimits)
           --retry-connrefused         retry even if connection is refused
           --retry-on-http-error=ERRORS    comma-separated list of HTTP errors to retry
      -O,  --output-document=FILE      write documents to FILE
      -nc, --no-clobber                skip downloads that would download to
                                         existing files (overwriting them)
    
      -c,  --continue                  resume getting a partially-downloaded file
           --progress=TYPE             select progress gauge type
           --show-progress             display the progress bar in any verbosity mode
      -N,  --timestamping              don't re-retrieve files unless newer than
                                         local
           --no-if-modified-since      don't use conditional if-modified-since get
                                         requests in timestamping mode
           --no-use-server-timestamps  don't set the local file's timestamp by
                                         the one on the server
       -T,  --timeout=SECONDS           set all timeout values to SECONDS
           --dns-timeout=SECS          set the DNS lookup timeout to SECS
           --connect-timeout=SECS      set the connect timeout to SECS
           --read-timeout=SECS         set the read timeout to SECS
      -w,  --wait=SECONDS              wait SECONDS between retrievals
           --waitretry=SECONDS         wait 1..SECONDS between retries of a retrieval
           --random-wait               wait from 0.5*WAIT...1.5*WAIT secs between retrievals
    
           --limit-rate=RATE           limit download rate e.g. 1M  1 MB/s      
    """
    pass
    logging.debug(f'wget id {self.runid}')
    loglev = LOGLEVELS[self.log.getEffectiveLevel()]
    cmd = ['wget',
           '--limit-rate', rate,
           '--continue', 
           '-O', f'{destpath}',
           f'{srcurl}']
    cmdstr = " ".join(cmd)
    logging.debug(f"wget command: {cmdstr} running...")
    cp = subprocess.run(cmd)
    self.log.debug(
        f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == "0":
        return self.runid
    else:
        self.log.error(f'non-zero return code for runid {self.runid}')
        return None
    


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
    
