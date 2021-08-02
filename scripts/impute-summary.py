#!/usr/bin/env python
#
#  Consume impute.tsv DF and emit summary numbers:
#  How many technologies? 
#  How many projects for each technology
#  How many runs for each project?
#

import argparse
import logging
import os
import sys

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import *

TECH_VERSIONS=['smartseq','10xv1','10xv2','10xv3']

def print_dict(title, d):
    print(title)
    for k in d.keys():
        print(f'{k}\t{d[k]}')
    print('')


def summarize(imputepath):
    idf = load_df(imputepath)
    logging.info(idf)
    tech_projnum = {}
    tech_projlists = {}
      
    for tech in TECH_VERSIONS:
        try:
            projlist  = list(idf.groupby('tech_version').proj_id.unique()[tech])
            numproj = len(projlist)
            tech_projnum[tech] = numproj
            
            for proj_id in projlist:
                pass
                
            
            
        except KeyError:
            pass
    print_dict('TECHNOLOGY PROJECT TOTALS',tech_projnum)
    
    


if __name__ == "__main__":

    # CLI setup
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
    parser.add_argument('imputetsv' )
    args = parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    ################## 
    
    imputefile = os.path.expanduser(args.imputetsv)
    logging.debug(f'consuming impute file: {imputefile}')
    out = summarize(imputefile)
    
    