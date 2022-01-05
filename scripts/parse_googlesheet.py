#!/usr/bin/env python
#
#  
#
#
#
import io
import os
import sys
import json
import pandas as pd
import requests
from urllib import parse
import logging

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *
from scqc.common import *
from scqc.esearch import *

def get_svensson_data():
    # google sheet URL 

    # as of 1/4/2022 this gives 117 experiments
    # 102 GEO ids
    # 1 SRP id
    # 4 PRJ ids
    # 5 SCP ids (single cell portal - broad inst)
    # 1 Nemo
    # 4 others (DOI/figshare/SCR/syn)
    url = "https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/export?gid=0&format=tsv"

    r = requests.post(url) 
    with io.StringIO(r.text) as imf:
        df = pd.read_csv(imf, sep='\t')

    # filter dataset for species
    ind = df.Organism.str.contains('Mouse')
    ind[ind.isna()] = False
    df = df.loc[ind,:]

    # filter for RNA seq
    ind = df.Measurement.str.contains('RNA-seq')
    ind[ind.isna()] = False
    df = df.loc[ind,:]

    # filter for brain
    ind = df.Tissue.str.contains('Brain')
    ind[ind.isna()] = False
    df = df.loc[ind,:]

    # make sure theres a data source avaiable
    ind = ~ df['Data location'].isna()
    # ind[ind.isna()] = False
    df = df.loc[ind,:]
    
    return df.reset_index(drop=True)

### Updated esearch query - allows for geo query 
def build_url(db = 'geo',
        species=["mus musculus"], 
        strategy=["rna seq"], 
        textword=["single cell","brain"],
        geo_acc =['GSE70844']
    ):
    """
    Given 3 lists, build search string for URL with following logic:
    
    (sp1 OR sp2) AND (st1 or st2) AND (tw1 AND tw2) 

    """
    baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db='

    DB_TAG ={
        'geo' : 'gds',
        'sra' : 'sra'
    }
    # GSE70844
    SRATAGS = { 'species': '[Organism]',
            'strategy' : '[Strategy]',
            'textword' : '[Text Word]'
            }
    GEOTAGS ={
            'geo_acc' : ['GEO Accession'],
            'species' : ['Organsim']
            }

    searchroot = f'{baseurl}{DB_TAG[db]}'
    tstrat = []
    tspec = []
    ttxt = []
    tgeo =[]
    

    if db =='sra':
        TAGS = SRATAGS
        for r in strategy:
            tstrat.append(f"{r}{TAGS['strategy']}")
        for t in textword:
            ttxt.append(f"{t}{TAGS['textword']}")
        tgeo =[]
        
    elif db =='geo':
        TAGS =GEOTAGS
        for g in geo_acc:
            tgeo.append(f"{g}{TAGS['geo_acc']}")
        tstrat =[]
        ttxt=[]

    for s in species:
        tspec.append(f"{s}{TAGS['species']}")


    spc = ' OR '.join(tspec)
    strat = ' OR '.join(tstrat)
    txtwrd = ' AND '.join(ttxt)
    geoid = ' OR '.join(tgeo)

    params = f"({spc}) AND ({strat}) AND ({txtwrd}) AND {geoid}"
    # self.log.debug(≥'uncoded params = {params}')
    encoded = parse.quote_plus(params)        
    furl = f"{searchroot}&term={encoded}"

    return furl    


def main():
    df = get_svensson_data()
    datalocs = df['Data location']
    geo_ids = datalocs.loc[datalocs.str.contains('GSE')].reset_index(drop=True)
    sra_ids = datalocs.loc[datalocs.str.contains('RP')].reset_index(drop=True)
    prj_ids = datalocs.loc[datalocs.str.contains('PRJ')].reset_index(drop=True)
    scp_ids = datalocs.loc[datalocs.str.contains('SCP')].reset_index(drop=True)

    t = [ i.split(',') for i in list(geo_ids)]    
    geo_ids = [item.strip() for sublist in t for item in sublist]


    


    all_sra_furl = build_url(db = 'sra',
        species=["mus musculus"], 
        strategy=["rna seq"], 
        textword=["single cell","brain"]
    )

    s_all = SraSearch(get_default_config())
    geo_uids = query_GEO()

    



def query_GEO(geo_ids):

    geo_furl = build_url(db = 'geo',
        species=["mus musculus"], 
        strategy=["rna seq"], 
        textword=["single cell","brain"],
        geo_acc =geo_ids)
    
    query_max = 50000
    query_start =0 
    all_uids = []
    while True:
        try :
            rgeo = requests.get(f'{geo_furl}&retstart={query_start}&retmax={query_max}&retmode=json')
            ergeo = json.loads(rgeo.content.decode('utf-8'))
            idlist = ergeo['esearchresult']['idlist']
            if len(idlist) >0 :
                for id in idlist:
                    all_uids.append(id)
                query_start += query_max
            else :
                break
        except Exception as ex:
            print(ex)
            break
            # log.error(traceback.format_exc(None))


    fetch_geo(all_uids)

    pass



def query_SRA(): # also works for PRJ
    pass

def query_SCP():
    # highly doubt we can get directly from SCP...
    pass

def fetch_geo(geo_uids) :
    batchsize = 100
    retstart = 0
    donelist = []
    
    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds'
    full_url = f"{efetch_url}&id={','.join(geo_uids[retstart:(retstart+batchsize)])}&retmode=doc"

    r =  requests.post(full_url)
    rd = r.content.decode().splitlines()
    
    # root = et.fromstring(rd)


    # SRA Run Selector -> link to traces -> extract prj id    

#print(df)
# locations = list(df['Data location'][df['Data location'].notna()])
# locations.sort()
# for s in locations:
#     print(f"{s}")