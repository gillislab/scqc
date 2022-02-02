#!/usr/bin/env python
import os
import sys
import pandas as pd
import logging


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.sra import Query
from scqc.utils import *

logging.getLogger().setLevel(logging.DEBUG)

# Query all 
q = Query(get_default_config())
q.tempdir = '/data/johlee/cross_mammal_xci/sra'
q.cachedir = '/data/johlee/cross_mammal_xci/cache'
q.metadir = '/data/johlee/cross_mammal_xci/metadata'

rdf = load_df(f'{q.metadir}/runs.tsv')
species_list = ['Bos taurus','Ovis aries','Rattus norvegicus','Macaca mulatta','Sus scrofa']
for spec in species_list:
    spec2 = spec.replace(' ','_')
    destdir= f'{q.tempdir}/sra/{spec2}'
    try:
        os.makedirs(destdir)
    except FileExistsError:
        pass

    srcurls = rdf.loc[rdf.sciname == spec ,['run_id','file_url']]
    for i in range(srcurls.shape[0]):
        try:
            download_wget(srcurls.file_url[i], destpath=f'{destdir}/{srcurls.run_id[i]}.sra',rate = '50M')
        except KeyError:
            pass