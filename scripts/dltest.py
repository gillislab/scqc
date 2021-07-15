#!/usr/bin/env python
#
# Test script for testing download methods. 
#
#  https://sra-download.ncbi.nlm.nih.gov/traces/era22/ERR/ERR4822/ERR4822769
#  https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR4822769/ERR4822769


import argparse
import requests
import logging
import time
import traceback


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
