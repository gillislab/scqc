
import os 
import glob
import sys
import subprocess
import argparse

import tarfile
import pandas as pd

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)

from scqc.utils import *
from scqc.star import *


FASTQS = pd.read_csv('/home/johlee/git/scqc/resource/biccn_fastq_paths.txt',header=None)
FASTQS.columns=['srcurl']
FASTQS = FASTQS.loc[ ~ FASTQS.duplicated(), :]  # just in case

CONF = get_default_config()

class AlignBICCN(object) :

    def __init__(self, config):
        self.config = config
        self.ar = AlignReads(config)
        self.clean  = self.config.get('star' ,'nocleanup')


    def execute(self,run_id ) : 
           
        # needs to be able to deal with multiple lanes. 
        tarpaths =  glob.glob(f'/data/biccn/nemo/{run_id}*fastq.tar')
        
        # unpack tars  # TODO multi thread
        # occasional issue?? NoneType?
        for tarpath in tarpaths :
            if os.path.isfile(tarpath) :
                print('untarring files' )
                self.untar_file(tarpath)

        # list the fastqs found
        
        subdirs = [ tarpath.replace('.fastq.tar','') for tarpath in tarpaths ] 
        fastqs = [ glob.glob(f'{subdir}/*.fastq.gz') for subdir in subdirs]
        fastqs = list(itertools.chain.from_iterable(fastqs))
        
        # What technology is the url suggesting?
        srcurls = [ f for f in FASTQS.srcurl  if run_id in f]
        is10xv2 = ['10X_V2'  in srcurl.upper()  for srcurl in srcurls ] 
        is10xv3 = ['10X_V3'  in srcurl.upper()  for srcurl in srcurls ] 
        isSSv4 = ['SSv4'  in srcurl.upper()  for srcurl in srcurls ] 

        if all(is10xv2) : 
            # which fastqs are the cDNA/barcodes?
            print(f'starting 10xv2 processing for {run_id}')
            barcodes = []
            cDNAs = []
            for fq in fastqs :
                l = self.get_read_length(fq)
                if l == 26 : 
                    barcodes.append(fq)
                if l > 30 :
                    cDNAs.append(fq)
            barcodes = ','.join(barcodes)
            cDNAs = ','.join(cDNAs)


            outfile_prefix = self.ar._run_star_10x(run_id = run_id, tech = '10xv2' , 
                bio_readpath = cDNAs, tech_readpath=barcodes , gzipped=True)
            if self.clean : 
                for fq in fastqs: 
                    os.remove(fq) 
            print(outfile_prefix)

        elif  all(is10xv3) : 
            print(f'starting 10xv3 processing for {run_id}')
            # which fastqs are the cDNA/barcodes? 
            barcodes = []
            cDNAs = []
            for fq in fastqs :
                l = self.get_read_length(fq)
                if l == 28 : 
                    barcodes.append(fq)
                if l > 30 :
                    cDNAs.append(fq)
            barcodes = ','.join(barcodes)
            cDNAs = ','.join(cDNAs)

            outfile_prefix = a._run_star_10x(run_id = run_id, tech = '10xv2' , 
                bio_readpath = cDNAs, tech_readpath=barcodes , gzipped=True)
             
        elif 'SSv4'  in srcurl.upper() : 
            # continue as smartseq. 
            #make manifest...
            pass


    def untar_file(self, tarpath):
        # ex /data/biccn/nemo/L8TX_180221_01_B10.fastq.tar
        # with tarfile.open(tarpath) as f : 
        #     f.extractall(os.path.dirname(tarpath) )
        cmd = ['tar', '-C',os.path.dirname(tarpath), '-xvf', tarpath ]
        run_command(cmd)


        if self.clean:
            os.remove(tarpath)

        # return(stderr, stdout, returncode)

    def get_read_length(self, gzpath):
        if os.path.isfile(gzpath):
            z_head = subprocess.Popen( f'zcat {gzpath} | head',
                stdout=subprocess.PIPE, shell=True ).communicate()
            
            stdout = z_head[0].decode().split('\n')
            l = len(stdout[1])

            return(l)

        
    
if __name__ == '__main__' :
    
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

    parser.add_argument('-c', '--config',
                            action="store",
                            dest='conffile',
                            default='~/git/scqc/etc/scqc.conf',
                            help='Config file path [~/git/scqc/etc/scqc.conf]')

    parser.add_argument('-r', '--run_id',
                        metavar='run_id',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='unique identifier for the run')



    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)


    if args.conffile is not None:
        cp = ConfigParser()
        cp.read(os.path.expanduser(args.conffile)) 
    else:
        cp =get_default_config()

    if args.run_id is not None:
        
        for run_id in args.run_id:            
            ab = AlignBICCN(cp)
            ab.execute(run_id)
