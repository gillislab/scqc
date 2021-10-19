
# from typing import Mapping

import pandas as pd
import glob
import os
import sys
gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *

# read in the manifest files
# NOTE these do not include the 
# manipath = '/data/biccn/nemo/manifest_all.tsv'
# metadatapath = '/data/biccn/nemo/metadata_all.tsv' 
# manifest = pd.read_csv(manipath, sep="\t")
# metadata = pd.read_csv(metadatapath, sep="\t")

# def main(manifest,metadata ):

#     # sdf = metadata 
    
#     urls = manifest.urls
#     urls = urls.str.split(',')

#     df = pd.DataFrame([url[0].split('/') for url in urls ])
#     df = df.iloc[:, 5:]
#     df.columns  = ['grant','author','modality','subspecimen_type','technique','species','data_type','region', 'file']
    

#     # extract the meta data from the url
# http_urls = mf.url


# def make_sdf(manifest, metadata ) :
#     # project data frames
#     metadata.project_id 



def make_dfs(dirname='/data/biccn/nemo/sc_mm_brain_manifests-20210923'):
    datafiles = glob.glob(f'{dirname}/*.tsv')
    metadatapaths = [ f for f in datafiles if '_metadata_' in f ] 
    manifestpaths = [ f for f in datafiles if f not in metadatapaths ]

    allmd = pd.DataFrame()
    for f in metadatapaths:
        tmpdf = pd.read_csv(f,sep='\t')
        allmd = pd.concat([allmd,tmpdf],axis=0)

    allmd = allmd.drop_duplicates().reset_index(drop=True)

    allattr = []
    attr_keys = ['sample_donor_id','sample_anatomical_region','sample_subspecimen_type','project_lab','project_organization']
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
    sdf['proj_id'] = allmd.project_id
    sdf['submission_id'] = allmd.project_grant

    edf =pd.DataFrame()
    edf['exp_id'] = allmd.sample_id
    edf['ext_ids']	=''
    edf['strategy']	='RNA-Seq'
    edf['source'] ='TRANSCRIPTOMIC'
    edf['lcp'] = allmd.study_full_name
    edf['samp_id'] = allmd.sample_id
    edf['proj_id'] = allmd.project_id
    edf['submission_id'] = allmd.project_grant

    
    tmpdf = allmd[['project_id','project_grant']].drop_duplicates().reset_index(drop=True)
    pdf = pd.DataFrame()
    pdf['proj_id'] = tmpdf.project_id
    pdf['ext_ids'] = ''
    pdf['title'] = ''
    pdf['abstract'] = ''
    pdf['submission_id'] = tmpdf.project_grant

    manifest = pd.DataFrame()
    for f in manifestpaths:
        tmpdf = pd.read_csv(f,sep='\t')
        manifest = pd.concat([manifest,tmpdf],axis=0)
    manifest[ ~ manifest.md5.isna() ]
    manifest=manifest.drop_duplicates().reset_index(drop=True)

    rdf2sdf = pd.merge(manifest[['file_name','size','sample_id']], sdf ,left_on = 'sample_id',right_on ='samp_id',how ='left')

    rdf = pd.DataFrame()
    rdf['run_id']  =  (rdf2sdf.file_name+rdf2sdf.file_name).replace('tar','')
    rdf['ext_ids']=''
    rdf['tot_spots']=''
    rdf['tot_bases']=''
    rdf['run_size']= ''
    rdf['publish_date']=''
    rdf['taxon']='10090'
    rdf['organism']= 'Mus musculus'
    rdf['nreads']=''
    rdf['basecounts']= ''	
    rdf['file_url']=[a[1] for a in manifest.urls.str.split(',')]
    rdf['file_size']=rdf2sdf.size
    rdf['exp_id'] =rdf2sdf.sample_id
    rdf['samp_id']=rdf2sdf.sample_id
    rdf['proj_id']=rdf2sdf.proj_id
    rdf['submission_id'] =rdf2sdf.submission_id

    return(rdf, edf, sdf, pdf)
    

rdf, edf, sdf, pdf =  make_dfs()

merge_write_df(rdf,'/data/hover/scqc/metadata/biccn_runs.tsv')
merge_write_df(edf,'/data/hover/scqc/metadata/biccn_experiments.tsv')
merge_write_df(sdf,'/data/hover/scqc/metadata/biccn_samples.tsv')
merge_write_df(pdf,'/data/hover/scqc/metadata/biccn_projects.tsv')