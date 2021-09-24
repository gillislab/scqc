#!/usr/bin/env python
#
# download nemo data from manifests... 
#
# file_name    md5    size    filetype    urls    sample_id    component_files

import argparse
import requests
import logging
import os
import sys
import time
import traceback

import pandas as pd
import numpy as np
import traceback


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import download_wget

def parse_urllist(urls, prot='http'):
    out = None
    for f in urls.split(','):
        if f.startswith(prot):
            out = f
    logging.debug(f'found url with prot={prot}: {f}')
    return out



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

    parser.add_argument('-m', '--manifest', 
                        type=str,
                        required=True, 
                        help='a NeMO manifest file.')   

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    df = pd.read_csv(args.manifest, index_col=0, sep='\t', comment="#")
    logging.info(f'got list of {len(df)} file_urls.')
    
    for index, row in df.iterrows():
        md5 = row['md5']
        fname=row['file_name']
        fsize = int(row['size'])
        urls = row['urls']
        logging.debug(f"fname={fname} md5={md5} fsize={fsize} urls={urls}")
        url = parse_urllist(urls)
        logging.debug(url)
        download_wget(url, f'./{fname}')
        time.sleep(2)
    