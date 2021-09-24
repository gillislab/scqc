#!/usr/bin/env python
#
#  Takes NeMO manifest.tsv file and outputs format usable by portal_client
# Includes only https URLS. 
#
#  MANIFEST
#      component_files    filetype    md5    sample_id    size    urls
# 0    gs://nemo-public/biccn-unbundled/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOs/SM-D9EPH_S23_E1-50_R1.fastq.gz;gs://nemo-public/biccn-unbundled/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOs/SM-D9EPH_S23_E1-50_R2.fastq.gz;NA;NA;NA;NA    FASTQ    3def8f8f18430ce38e5f377dc668ebb1    553517840    159836160    gs://nemo-public/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOs/SM-D9EPH_S23_E1-50.fastq.tar,https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/MOs/SM-D9EPH_S23_E1-50.fastq.tar
#
# PORTAL CLIENT 
# id    md5    size    urls
#  8d65351c0cb6f7eb9726a1c463f1c34e    626b3b87b8958e1db84489d727d16607    2353731    http://downloads.hmpdacc.org/data/HM16STR/HMDEMO/SRP002429/stool/not_affected/SRS066784.fsa
#
#
#

import os
import sys
#gitpath=os.path.expanduser("~/git/cshlwork")
#sys.path.append(gitpath)

import pandas as pd
import numpy as np
import traceback

if not len(sys.argv) > 1:
    print("need input file arg(s)")
    sys.exit()

newdf = None  
  
infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    # print(f"{infile}")
    try:
        df = pd.read_csv(infile, index_col=0, sep='\t', comment="#")
        newdf = df[['md5', 'md5', 'size','urls']]
        newdf.drop_duplicates(inplace=True)
        newdf = newdf.reset_index(drop=True)
        newdf.to_csv(sys.stdout, sep='\t')     
   
    except Exception as ex:
        print(traceback.format_exc(None))
        raise ex

