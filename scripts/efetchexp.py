#!/usr/bin/env python
#
# Prints experiment XML data 
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=ERX4692753
#

import argparse
import requests
import logging
import time
import traceback

import xml.etree.ElementTree as et

from bs4 import BeautifulSoup as bs


sra_efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra'
query_sleep = 0.5

def query_experiment_package_set( xid):
    """
    Query XML data for this experiment ID. 

    """
    log = logging.getLogger()
    xmldata = None
    try:
        url = f"{sra_efetch}&id={xid}"
        log.debug(f"fetch url={url}")

        while True:
            r = requests.post(url)
            if r.status_code == 200:
                xmldata = r.content.decode()
                log.debug(f'good HTTP response for {xid}')
                break
            else:
                log.warn(
                    f'bad HTTP response for id {xid}. retry in 10s')
                time.sleep(10)

    except Exception as ex:
        log.error(f'problem with NCBI id {xid}')
        logging.error(traceback.format_exc(None))

    finally:
        log.debug(
            f"sleeping {query_sleep} secs between fetch calls...")
        time.sleep(query_sleep)
    return xmldata
    
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
    
    parser.add_argument('exp_ids',
                        nargs='+'
                        )
    
    args = parser.parse_args()
    
    for exp in args.exp_ids:
        xmlstr = query_experiment_package_set(exp)
        #root = et.fromstring(xmlstr)
        #print(et.tostring(root, encoding='unicode', pretty_print=True))
        print(bs(xmlstr, "xml").prettify())
        