
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



def make_dfs(dirname='/data/biccn/sc_mm_brain_manifests-20211019'):
    datafiles = glob.glob(f'{dirname}/*.tsv')
    metadatapaths = [ f for f in datafiles if '_metadata_' in f ] 
    manifestpaths = [ f for f in datafiles if f not in metadatapaths ]

    allmd = pd.DataFrame()
    for f in metadatapaths:
        tmpdf = pd.read_csv(f,sep='\t')
        allmd = pd.concat([allmd,tmpdf],axis=0)

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

    manifest = pd.DataFrame()
    for f in manifestpaths:
        tmpdf = pd.read_csv(f,sep='\t')
        manifest = pd.concat([manifest,tmpdf],axis=0)
    manifest = manifest[ ~ manifest.md5.isna() ]
    manifest=manifest.drop_duplicates().reset_index(drop=True)

    rdf2sdf = pd.merge(manifest[['file_name','size','sample_id','urls']], sdf ,left_on = 'sample_id',right_on ='samp_id',how ='left')

    rdf = pd.DataFrame()
    rdf['run_id']  =  (rdf2sdf.proj_id+'_'+rdf2sdf.file_name).str.replace('.fastq.tar','',regex=False)
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

    return(rdf, edf, sdf, pdf)
    

rdf, edf, sdf, pdf =  make_dfs()
print(rdf.shape, edf.shape, sdf.shape, pdf.shape)
outdir = '/data/johlee'
# outdir = '/data/hover/scqc/metadata'
merge_write_df(rdf,f'{outdir}/biccn_runs.tsv')
merge_write_df(edf,f'{outdir}/biccn_experiments.tsv')
merge_write_df(sdf,f'{outdir}/biccn_samples.tsv')
merge_write_df(pdf,f'{outdir}/biccn_projects.tsv')