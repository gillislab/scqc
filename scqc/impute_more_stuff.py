
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
def get_vdb_data(runfile, outfile):
    runs = pd.read_csv(runfile, sep="\t")

    allrows =list()
    for run_id in runs.run_id : 
        print(run_id)
        cmd = ['vdb-dump', 
            '--rows', '1',
            '--columns', 'READ_LEN,SPOT_COUNT',
            run_id]
        stderr, stdout, rc = run_command(cmd)
        # print("...",stdout.strip())
        tmp = stdout.strip().split('\n')
        row = [i.split(': ' )[1].strip() for i in   stdout.strip().split('\n') ]
        row.append(run_id)

        allrows.append(row)

    run_read_lengths = pd.DataFrame(allrows,columns= ['read_length','nspots','runid'])
    run_read_lengths.to_csv(outfile, sep="\t")


if __name__== "__main__":

    parser = argparse.ArgumentParser()


    parser.add_argument('-r','--runfile',
                        metavar='runfile',
                        type=str,
                        default=None,
                        help='run files')

    parser.add_argument('-o', '--outfile',
                        metavar='outfile',
                        type=str,
                        default=None,
                        help='Outfile. ')

    args = parser.parse_args()

    get_vdb_data(args.runfile, args.outfile)

    