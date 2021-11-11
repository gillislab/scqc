#!/usr/bin/env python
import os
import sys
import ast
import glob

import numpy as np
import scanpy as sc
import pandas as pd



gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *


h5paths = glob.glob('/data/hover/scqc/output/*.h5ad')

obs_scalar_columns = ['n_genes_by_counts','total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_500_genes',
    'pct_counts_mt','pct_counts_ribo','pct_counts_female',  'pct_counts_male',  
    'pct_counts_essential','pct_counts_highly_variable','pct_counts_housekeeping','gini','corr_to_mean','max_corr_with_others']

obs_discrete_columns = ['tech','class_label','subclass_label']
alldf = pd.DataFrame()
allsampattr =pd.DataFrame()

for h5path in h5paths:
    adata = sc.read_h5ad(h5path)
    
    proj_id = os.path.basename(h5path).replace('.h5ad','')
    print(proj_id)
    title = adata.uns['title']
    ncells = adata.shape[0]
    
    obs_dat = adata.obs[obs_scalar_columns].mean(0)
    prop_female = obs_dat.pct_counts_female / (obs_dat.pct_counts_female+obs_dat.pct_counts_male)
    prop_male = obs_dat.pct_counts_male / (obs_dat.pct_counts_female+obs_dat.pct_counts_male)


    # sample_attributes counts
    tmp = pd.merge(adata.obs[['cell_id','samp_id']], adata.uns['samp_df'],on = 'samp_id',how='left')
    asdict = tmp['attributes'].reset_index(drop=True).apply(ast.literal_eval)

    try: 
        source_names_cell_counts = asdict.apply( lambda x :x['source_name'] ).value_counts()
        source_names_cell_counts.index=['source_name: ' +s  for s in list(source_names_cell_counts.index)]
    except KeyError:
        print(f'no source name for {proj_id}')



    df = pd.DataFrame({proj_id:obs_dat})
    df = df.append(pd.DataFrame({ proj_id: [prop_female, prop_male] }, index=['prop_female','prop_male'] ) )

    # append discrete counts
    df = df.append(  pd.DataFrame({proj_id:  adata.obs[obs_discrete_columns].melt().value_counts()  } ))

    #append source_name counts
    df = df.append(  pd.DataFrame({proj_id: source_names_cell_counts   } ))


    # for cleaning
    alldf = pd.merge(alldf,df,right_index=True,left_index=True, how='outer')
    allsampattr = allsampattr.append(pd.DataFrame( adata.uns['samp_df'].attributes))

    print(alldf)


# all source names 
# AON -> anterior olfactory neuron
# V-SVZ -> ventricular-subventricular zone 

# ['source_name: AON viral injection',
#  'source_name: Adult postnatal V-SVZ tissue',
#  'source_name: Aged brain endothelial cells',
#  'source_name: Blood',
#  'source_name: Bone flap',
#  'source_name: Brain cortex',
#  'source_name: Brain meninges',
#  'source_name: Brain neuron',
#  'source_name: Brain tissue',
#  'source_name: Brain',
#  'source_name: CD11+ myeloid cells from tumor-bearing cerebellum treated with GDC-0449',
#  'source_name: CD11+ myeloid cells from tumor-bearing cerebellum treated with radiation',
#  'source_name: CD11+ myeloid cells from tumor-bearing cerebellum',
#  'source_name: CD45+ cells from tumor-bearing cerebellum treated with radiation',
#  'source_name: CD45+ cells',
#  'source_name: CD45+ dura single cells',
#  'source_name: CD45-+ brain single cells',
#  'source_name: Cerebellum',
#  'source_name: Control: 10pg Takara Control Total RNA',
#  'source_name: Corpus Callosum',
#  'source_name: Corpus_Callosum',
#  'source_name: Cortex',
#  'source_name: Cortical Astrocytes',
#  'source_name: Cuprizone_demyelination_corpus_callosum',
#  'source_name: Cuprizone_remyelination_corpus_callosum',
#  'source_name: Dentate Gyrus',
#  'source_name: E15.brain.cell',
#  'source_name: Early postnatal V-SVZ tissue',
#  'source_name: Embryonic V-SVZ tissue',
#  'source_name: Embryonic cortex and ganglionic eminence',
#  'source_name: Facial_nerve_axotomy_d14_Facial_nucleus',
#  'source_name: Facial_nerve_axotomy_d3_Facial_nucleus',
#  'source_name: Facial_nucleus',
#  'source_name: Forebrain, fresh tissue',
#  'source_name: Full Brain',
#  'source_name: Galea',
#  'source_name: Hindbrain, fresh tissue',
#  'source_name: Hippocampus',
#  'source_name: Immune cells',
#  'source_name: LPS_24h_Female',
#  'source_name: LPS_24h_Male',
#  'source_name: Late postnatal V-SVZ tissue',
#  'source_name: MJ5_Emx60ctx_KO_S1_I3',
#  'source_name: MJ5_Emx64ctx_Ctl_S2_I3',
#  'source_name: MJ5_Emx65ctx_Ctl_S3_I3',
#  'source_name: MJ5_Emx73ctx_KO_S4_I3',
#  'source_name: Micro dissected ventricular-subventricular zone',
#  'source_name: Mouse brain meninges',
#  'source_name: Murine brain tumor infiltrated cells',
#  'source_name: Olfactory bulb',
#  'source_name: P21 Brain',
#  'source_name: Pdgfra-Cre-LoxP-GFP',
#  'source_name: Pdgfra-H2B-GFP',
#  'source_name: Pons, fresh tissue',
#  'source_name: Pooled Cd11+ myeloid cells from 3 mouse cerebella',
#  'source_name: Pooled Cd11+ myeloid cells from peripheral blood of 3 mice',
#  'source_name: Primary Visual Cortex (VISp)',
#  'source_name: Rhombomere 1 cells',
#  'source_name: Saline_24h_Female',
#  'source_name: Saline_24h_Male',
#  'source_name: Single E15.5 FB mouse Th-eGFP+ sorted dopaminergic neuron',
#  'source_name: Single E15.5 MB mouse Th-eGFP+ sorted dopaminergic neuron',
#  'source_name: Single P7 FB mouse Th-eGFP+ sorted dopaminergic neuron',
#  'source_name: Single P7 MB mouse Th-eGFP+ sorted dopaminergic neuron',
#  'source_name: Single P7 OB mouse Th-eGFP+ sorted dopaminergic neuron',
#  'source_name: Spinal_Cord',
#  'source_name: TNFtg',
#  'source_name: WT_cerebellum_03_w',
#  'source_name: WT_cerebellum_16_w',
#  'source_name: WT_cerebellum_embryonal',
#  'source_name: WT_corpus_callosum_03_w',
#  'source_name: WT_corpus_callosum_16_w',
#  'source_name: WT_cortex_03_w',
#  'source_name: WT_cortex_16_w',
#  'source_name: WT_forebrain_embryonal',
#  'source_name: WT_hippocampus_03_w',
#  'source_name: WT_hippocampus_16_w',
#  'source_name: WT_midbrain_embryonal',
#  'source_name: WT_spinal_cord_03_w',
#  'source_name: WT_spinal_cord_16_w',
#  'source_name: WT_spinal_cord_embryonal',
#  'source_name: WT',
#  'source_name: Whole brain',
#  'source_name: adult whole kidney',
#  'source_name: brain immune cells FACS sorted CD45+',
#  'source_name: brain slice microglia',
#  'source_name: brain tumor xenograft',
#  'source_name: brain-hypothalamus',
#  'source_name: brain; tumor',
#  'source_name: brain',
#  'source_name: freshly isolated brain endothelial cells',
#  'source_name: hippocampal cells',
#  'source_name: in vitro cultured microglia',
#  'source_name: lateral hypothalamic area',
#  'source_name: meninges',
#  'source_name: microglia',
#  'source_name: mouse brain cortex',
#  'source_name: photoconverted CD4+ T cells isolated from inguinal lymph node + mesenteric lymph node + spleen + brain',
#  'source_name: photoconverted CD4+ T cells isolated from inguinal lymph node + spleen + brain',
#  'source_name: photoconverted CD4+ T cells isolated from mesenteric lymph node + spleen + brain',
#  'source_name: posterior Piriform cortex viral injection',
#  'source_name: single brain cell',
#  'source_name: striatal Rbpj-K -/- astrocytes',
#  'source_name: whole brain without cerebellum']