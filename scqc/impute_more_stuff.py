
import os 
import glob
import sys
import subprocess
import argparse

import pandas as pd


gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *






# look for 10x runs  
def get_vdb_data(start ,stop, outfile):

    print(f'starting - {start}:{stop} and saving at {outfile}')
    runs = glob.glob('/data/hover/scqc/cache2/sra/*.sra')
    runs = [run_id.replace('/data/hover/scqc/cache2/sra/','') for run_id in runs ]
    runs = [run_id.replace('.sra','') for run_id in runs ]

    allrows =list()
    for run_id in runs[start:stop] : 
        # print(run_id)
        cmd = ['vdb-dump', 
            '--rows', '1',
            '--columns', 'READ_LEN,SPOT_COUNT',
            run_id]
        stderr, stdout = run_command(cmd)
        # print("...",stdout.strip())
        tmp = stdout.strip().split('\n')
        row = [i.split(': ' )[1].strip() for i in   stdout.strip().split('\n') ]
        row.append(run_id)

        allrows.append(row)

    run_read_lengths = pd.DataFrame(allrows,columns= ['read_length','nspots','runid'])
    run_read_lengths.to_csv(outfile, sep="\t")


if __name__== "__main__":

    parser = argparse.ArgumentParser()


    parser.add_argument('-e','--stop',
                        metavar='stop',
                        type=int,
                        default=1,
                        help='stop index')

    parser.add_argument('-i','--start',
                        metavar='start',
                        type=int,
                        default=0,
                        help='start index')

    parser.add_argument('-o', '--outfile',
                        metavar='outfile',
                        type=str,
                        default=None,
                        help='Outfile. ')

    args = parser.parse_args()

    get_vdb_data(args.start,args.stop,args.outfile)
    