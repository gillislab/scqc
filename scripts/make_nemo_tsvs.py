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

def make_dfs(manifest, metadata):

    allmd = pd.read_csv(metadata,sep='\t')
    #allmd = pd.concat([allmd,tmpdf],axis=0)

    allmd = allmd.drop_duplicates().reset_index(drop=True)

    allattr = []
    attr_keys = ['sample_anatomical_region','sample_subspecimen_type','project_lab','project_organization'] #'sample_donor_id',
    for i in range(allmd.shape[0]):
        attr = {}
        for ky in attr_keys:
            attr[ky] = allmd.loc[i, ky]
        allattr.append(str(attr))

    sdf = pd.DataFrame()
    sdf['samp_id'] = allmd.sample_id
    sdf['ext_ids'] = ''
    sdf['taxon'] = '10090'
    sdf['sciname'] = 'Mus musculus'
    sdf['title'] = 	allmd.study_full_name
    sdf['attributes'] = allattr
    sdf['proj_id'] = allmd.project_id.str.replace(';','-')
    sdf['submission_id'] = allmd.project_grant

    edf =pd.DataFrame()
    edf['exp_id'] = allmd.sample_id
    edf['ext_ids']	=''
    edf['strategy']	='RNA-Seq'
    edf['source'] ='TRANSCRIPTOMIC'
    edf['lcp'] = allmd.study_full_name
    edf['samp_id'] = allmd.sample_id
    edf['proj_id'] = allmd.project_id.str.replace(';','-')
    edf['submission_id'] = allmd.project_grant

    
    tmpdf = allmd[['project_id','project_grant']].drop_duplicates().reset_index(drop=True)
    pdf = pd.DataFrame()
    pdf['proj_id'] = tmpdf.project_id.str.replace(';','-')
    pdf['ext_ids'] = ''
    pdf['title'] = ''
    pdf['abstract'] = ''
    pdf['submission_id'] = tmpdf.project_grant

    mdf = pd.read_csv(manifest,sep='\t')
    #manifest = pd.concat([manifest,tmpdf],axis=0)
    mdf = mdf[ ~ mdf.md5.isna() ]
    mdf=mdf.drop_duplicates().reset_index(drop=True)

    rdf2sdf = pd.merge(mdf[['file_id','size','sample_id','urls']], sdf ,left_on = 'sample_id',right_on ='samp_id',how ='left')

    rdf = pd.DataFrame()
    #rdf['run_id']  =  (rdf2sdf.proj_id+'_'+rdf2sdf.file_id).str.replace('.fastq.tar','',regex=False)
    rdf['run_id']  =  rdf2sdf.file_id.str.replace('.fastq.tar','',regex=False)
    rdf['ext_ids']=''
    rdf['tot_spots']=''
    rdf['tot_bases']=''
    rdf['run_size']= ''
    rdf['publish_date']=''
    rdf['taxon']='10090'
    rdf['organism']= 'Mus musculus'
    rdf['nreads']=''
    rdf['basecounts']= ''	
    rdf['file_url']=[a[1] for a in rdf2sdf.urls.str.split(',')]
    rdf['file_size']=rdf2sdf['size']
    rdf['exp_id'] =rdf2sdf.sample_id
    rdf['samp_id']=rdf2sdf.sample_id
    rdf['proj_id']=rdf2sdf.proj_id
    rdf['submission_id'] =rdf2sdf.submission_id

    # remove bulk samples that snuck in
    bulksamp = allmd.sample_id[allmd.sample_subspecimen_type=='Bulk']
    rdf = rdf.drop(rdf[rdf.samp_id.isin(bulksamp) ].index).reset_index(drop=True)
    sdf = sdf.drop(sdf[sdf.samp_id.isin(bulksamp) ].index).reset_index(drop=True)
    edf = edf.drop(edf[edf.samp_id.isin(bulksamp) ].index).reset_index(drop=True)

    # remove human samples that snuck in
    humansamp = allmd.sample_id[allmd.sample_organism=='Human']
    rdf = rdf.drop(rdf[rdf.samp_id.isin(humansamp) ].index).reset_index(drop=True)
    sdf = sdf.drop(sdf[sdf.samp_id.isin(humansamp) ].index).reset_index(drop=True)
    edf = edf.drop(edf[edf.samp_id.isin(humansamp) ].index).reset_index(drop=True)

    # remove projects that shouldn't be there
    validproj = set(rdf.proj_id)
    pdf = pdf[pdf.proj_id.isin(validproj)].reset_index(drop=True)
    
    # add data_source to all (nemo)
    rdf['data_source'] = 'nemo'
    sdf['data_source'] = 'nemo'
    edf['data_source'] = 'nemo'
    
    return(rdf, edf, sdf, pdf)




    
def dummy():
    rdf, edf, sdf, pdf =  make_dfs()
    print(rdf.shape, edf.shape, sdf.shape, pdf.shape)
    outdir = '/data/johlee'
    # outdir = '/data/hover/scqc/metadata'
    merge_write_df(rdf,f'{outdir}/nemo_runs.tsv')
    merge_write_df(edf,f'{outdir}/nemo_experiments.tsv')
    merge_write_df(sdf,f'{outdir}/nemo_samples.tsv')
    merge_write_df(pdf,f'{outdir}/nemo_projects.tsv')


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

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    rdf, edf, sdf, pdf =  make_dfs(manifest=args.manifest, 
                                   metadata=args.metadata)
    logging.info(f'got DFs: {pdf}\n{rdf}\n{edf}\n{sdf}')


