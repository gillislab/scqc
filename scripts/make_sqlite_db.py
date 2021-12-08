#!/usr/bin/env python
#
import argparse
import ast
import logging
import os
import sys
import traceback

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

import pandas as pd
import numpy as np
import traceback
import sqlite3 

from collections import defaultdict

from scqc.utils import *

IMAGE_PREFIX= "http://labshare.cshl.edu/shares/gillislab/resource/MetaQC/pcimage"

# Shiny URL form:
#  https://gillisweb.cshl.edu/MetaQC/?_inputs_&projid=%22SRP167086%22

DB_FIELDS = [   'proj_id', 
                'title',
                'submission_id', 
                'data_source', 
                'date',  
                'n_cells',  
                'n_smartseq',  
                'n_10xv2',  
                'n_10xv3',                 
                'gini_coeff', 
                'corr_to_mean',                  
                'class_pr',  
                'subclass_pr_hier', 
                'total_counts',  
                'ribo',  
                'mt',  
                'essential', 
                'cell_cycle', 
                'housekeeping',  
                'female',  
                'male',  
                'n_genes_by_counts',  
                'gini_cell',  
                'highly_variable',  
                'top_50_gene',  
                'top_100_gene', 
                'top_200_gene', 
                'top_500_gene' ]

TECH_FIELDS = ['smartseq', '10xv2','10xv3']

#
#  map from dataframe column to database column name
#
DF2DB = {
    'total_counts' : 'total_counts' , 
    'n_genes_by_counts' : 'n_genes_by_counts', 
    'corr_to_mean' : 'corr_to_mean', 
    'gini_x': 'gini_cell',
    'max_corr_with_others' : 'max_corr_with_others' , 
    'ribo':  'ribo', 
    'mt': 'mt' , 
    'essential':'essential' , 
    'cell_cycle':'cell_cycle' ,
    'housekeeping':'housekeeping' , 
    'female': 'female', 
    'male':'male' , 
    'highly_variable':'highly_variable' , 
    'top_50_gene' : 'top_50_gene',
    'top_100_gene' : 'top_100_gene', 
    'top_200_gene' : 'top_200_gene', 
    'top_500_gene' : 'top_500_gene', 
    'proj_id' : 'proj_id', 
    'ext_ids': 'ext_ids',
    'title' : 'title', 
    #'abstract' : 'abstract', 
    'submission_id' : 'submission_id', 
    'data_source' : 'data_source', 
    'class_pr' : 'class_pr',
    'subclass_pr': 'subclass_pr', 
    'subclass_pr_hier' : 'subclass_pr_hier', 
    'gini_y':'gini_coeff',  
    'n_cells': 'n_cells', 
#   'n_techs' :  'n_techs',
    'date' :  'date'      
    }

DB2DF = {
    'total_counts' : 'total_counts' , 
    'n_genes_by_counts' : 'n_genes_by_counts', 
    'corr_to_mean' : 'corr_to_mean', 
    'gini_cell': 'gini_x',
    'max_corr_with_others' : 'max_corr_with_others' , 
    'ribo':  'ribo', 
    'mt': 'mt' , 
    'essential':'essential' , 
    'cell_cycle':'cell_cycle' ,
    'housekeeping':'housekeeping' , 
    'female': 'female', 
    'male':'male' , 
    'highly_variable':'highly_variable' , 
    'top_50_gene' : 'top_50_gene',
    'top_100_gene' : 'top_100_gene', 
    'top_200_gene' : 'top_200_gene', 
    'top_500_gene' : 'top_500_gene', 
    'proj_id' : 'proj_id', 
    'ext_ids': 'ext_ids',
    'title' : 'title', 
    #'abstract' : 'abstract', 
    'submission_id' : 'submission_id', 
    'data_source' : 'data_source', 
    'class_pr' : 'class_pr',
    'subclass_pr': 'subclass_pr', 
    'subclass_pr_hier' : 'subclass_pr_hier', 
    'gini_coeff':'gini_y',  
    'n_cells': 'n_cells', 
#   'n_techs' :  'n_techs',
    'date' :  'date'      
    }




def def_value():
    return 0


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

    parser.add_argument('-i', '--infile', 
                        type=str,
                        required=True, 
                        help='A summary .tsv file. ')

    parser.add_argument('-b', '--dbfile', 
                        type=str,
                        required=True, 
                        help='Output SQLite DB file. ')    

    parser.add_argument('-t', '--tsvfile', 
                        type=str,
                        default=None,
                        required=False, 
                        help='Output .tsv file ')        

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    df = load_df(args.infile)
    logging.info(f'got list of {len(df)} datasets.')
    
    con = sqlite3.connect(args.dbfile)
    cur = con.cursor()
    #cur.execute('''CREATE TABLE datasets
    #(proj_id text, date text, n_cells integer, trans text, symbol text, qty real, price real)''')
    
    # density   double slider
    # % dropout double slider
    
    try:
        cur.execute('''create table datasets
                   ( proj_id text primary key,
                     title text,
                     submission_id text, 
                     data_source text, 
                     date text, 
                     n_cells integer, 
                     n_smartseq integer, 
                     n_10xv2 integer, 
                     n_10xv3 integer,                
                     gini_coeff real,
                     corr_to_mean real,                 
                     class_pr real, 
                     subclass_pr_hier real,
                     total_counts real, 
                     ribo real, 
                     mt real, 
                     essential real,
                     cell_cycle real,
                     housekeeping real, 
                     female real, 
                     male real, 
                     n_genes_by_counts real, 
                     gini_cell real, 
                     highly_variable real , 
                     top_50_gene real, 
                     top_100_gene real,
                     top_200_gene real,
                     top_500_gene real,
                     image text,
                     shinyurl text
                       )''')
        logging.info(f'table datasets created...')
        con.commit()

    except sqlite3.OperationalError as soe:
        logging.error(traceback.format_exc(None))
        logging.info('table exists?')
    
    
    tsvdata = []
            
    for index, row in df.iterrows():
        tl = []
        td = defaultdict(def_value)
        tech_dict = ast.literal_eval(row.n_techs)
        for k in tech_dict.keys():
            td[k] = tech_dict[k]
        for dbf in DB_FIELDS[:6]:
            tl.append(row[DB2DF[dbf]])
        
        # expand tech n_tech to 3 columns
        for tech in TECH_FIELDS:
            tl.append(td[tech])
        
        # continue with rest of DF columns
        for dbf in DB_FIELDS[9:]:
            tl.append(row[DB2DF[dbf]])
        
        # create image URL
        imageurl = f"{IMAGE_PREFIX}/{row.proj_id}_umap.png"
        tl.append(imageurl)
        
        # create shiny URL
        shinyurl = f'https://gillisweb.cshl.edu/MetaQC/?_inputs_&projid="{row.proj_id}"'
        tl.append(shinyurl)
        
        logging.debug(f'length = {len(tl)}') 
        tsvdata.append(tl)
        tt = tuple(tl)
        logging.debug(tt)
        
        try:
            inquery = f"insert into datasets values {tt}"
            logging.debug(inquery)
            cur.execute(inquery)
            con.commit()
        
        except Exception as e:
            logging.error(traceback.format_exc(None))

    con.close()        
    logging.info('done w/ dbfile.')
    
    if args.tsvfile is not None:
        tsvfields = DB_FIELDS.copy()
        tsvfields.append('image')
        tsvfields.append('shinyurl')
        df = pd.DataFrame(data=tsvdata, columns=tsvfields)
        logging.debug(f"dataframe: \n{df}")
        merge_write_df(df, args.tsvfile )
    
    
    