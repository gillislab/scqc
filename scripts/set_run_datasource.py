#!/usr/bin/env python
# from typing import Mapping
import argparse
import pandas as pd
import glob
import os
import sys
gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import *
from scqc.common import *
from scqc.nemo import *

#
# Set data_source based on value of project id, e.g.
#   rdf['data_source'] = np.where( rdf['proj_id'] == 'U19_Huang-Arlotta','nemo',rdf['data_source'])
#


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

    parser.add_argument('-f', '--manifest', 
                        metavar='manifest',
                        required=True, 
                        type=str, 
                        help='a (fixed) Nemo manifest TSV. ')
       
    parser.add_argument('-t', '--metadata', 
                        metavar='metadata',
                        required=True,  
                        type=str, 
                        help='a Nemo metadata file TSV. ') 
    
    parser.add_argument('-p', '--outprefix',
                        metavar='outprefix',
                        required=False,
                        type=str,
                        default='./',
                        help='prefix for all output dataframe TSVs'
                        )

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    make_dfs( manifest=args.manifest, 
              metadata=args.metadata, 
              prefix=args.outprefix)
