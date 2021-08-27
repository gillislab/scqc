#!/usr/bin/env python
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

BICCN_PROJS = pd.read_csv('/home/johlee/git/scqc/resource/biccn_projects.tsv',index_col=0,sep="\t")

CONF = get_default_config()

class AlignBICCN(object) :

    def __init__(self, config):
        self.config = config
        self.nocleanup  = self.config.get('star' ,'nocleanup') # if false, cleans tmp dir
        self.log = logging.getLogger('biccn_align')
        self.tempdir = os.path.expanduser(
            self.config.get('star', 'tempdir')) 



    def execute(self,run_id ) : 
           
        # needs to be able to deal with multiple lanes. 
        ar = AlignReads(self.config)

        tarpaths =  glob.glob(f'/data/biccn/nemo/{run_id}*fastq.tar')
        if len(tarpaths) ==0 :
            # no fastqs with that id found. check the list if it should be there. 
            rdf = BICCN_PROJS.loc[BICCN_PROJS.run_id ==run_id,:]
            if rdf.shape[0] == 0:
                raise NameError ('Invalid run_id')
            else :
                # get the files that we do have...
                fqs_avail = os.listdir('/data/biccn/nemo/')

                # get the ones that we expect to have
                fqs_required = rdf.uid
                tarpaths = [f'/data/biccn/nemo/{uid}.fastq.tar' for uid in rdf.uid if f'{uid}.fastq.tar' in fqs_avail ]
                if len(tarpaths) == 0 :
                    raise NameError (f'no runs downloaded yet for {run_id}')

                # Require all of the runs to be downloaded. If not, skip.
                # elif len(tarpaths) < len(rdf.uid) : 
                #     raise NameError (f'not all runs downloaded yet for {run_id}')
           
        # unpack tars  # TODO multi thread
        # occasional issue?? NoneType?
        for tarpath in tarpaths :
            if os.path.isfile(tarpath) :
                print('untarring files' )
                self.untar_file(tarpath)

        # list the fastqs found
        
        subdirs = [ tarpath.replace('.fastq.tar','') for tarpath in tarpaths ] 
        fastqs = [ glob.glob(f'{subdir}/*.fastq.gz') for subdir in subdirs]
        ulfastqs = list(itertools.chain.from_iterable(fastqs))
        
        # What technology is the url suggesting?
        isSSv4  = [True if 'SSv4' in f else False for f in rdf.tech ]
        is10xv2 = [True if '10x_v2' in f else False for f in rdf.tech ]
        is10xv3 = [True if '10x_v3' in f else False for f in rdf.tech ]

        if all(is10xv2) : 
            # which fastqs are the cDNA/barcodes?
            print(f'starting 10xv2 processing for {run_id}')
            barcodes = []
            cDNAs = []
            for fq in ulfastqs :
                l = self.get_read_length(fq)
                if l == 26 : 
                    barcodes.append(fq)
                if l > 30 :
                    cDNAs.append(fq)
            barcodes = ','.join(barcodes)
            cDNAs = ','.join(cDNAs)


            outfile_prefix = ar._run_star_10x(run_id = run_id, tech = '10xv2' , 
                bio_readpath = cDNAs, tech_readpath=barcodes , gzipped=True)
            if not self.nocleanup : 
                for fq in ulfastqs: 
                    os.remove(fq) 
            print(outfile_prefix)

        elif  all(is10xv3) : 
            print(f'starting 10xv3 processing for {run_id}')
            # which fastqs are the cDNA/barcodes? 
            barcodes = []
            cDNAs = []
            for fq in ulfastqs :
                l = self.get_read_length(fq)
                if l == 28 : 
                    barcodes.append(fq)
                if l > 30 :
                    cDNAs.append(fq)
            barcodes = ','.join(barcodes)
            cDNAs = ','.join(cDNAs)

            outfile_prefix = a._run_star_10x(run_id = run_id, tech = '10xv2' , 
                bio_readpath = cDNAs, tech_readpath=barcodes , gzipped=True)
            if not self.nocleanup : 
                for fq in ulfastqs: 
                    os.remove(fq) 
            print(outfile_prefix)

        elif all(isSSv4) : 
            # continue as smartseq. 
            # make manifest...
            manirows = [fqs if len(fqs) ==2 else fqs+['-'] for fqs in fastqs ] 
            manifest = pd.DataFrame(manirows ,columns =['read1','read2'])
            manifest['cell_id'] = subdirs
            manipath = f'{self.tempdir}/{run_id}_manifest.tsv'
            manifest.to_csv(manipath,sep="\t" ,header=None, index=None)

            outfile_prefix = ar._run_star_smartseq(run_id, manipath,gzipped= True)


            # if not self.nocleanup : 
            #     for fq in ulfastqs: 
            #         os.remove(fq) 
            print(outfile_prefix)


        else : 
            raise AssertionError (f'technology for the run is unclear for {run_id}')
            


    def untar_file(self, tarpath):
        # ex /data/biccn/nemo/L8TX_180221_01_B10.fastq.tar
        # with tarfile.open(tarpath) as f : 
        #     f.extractall(os.path.dirname(tarpath) )

        outdir = os.path.dirname(tarpath) 
        outbn = os.path.basename(tarpath).replace('.fastq.tar','')
        if os.path.isdir(f'{outdir}/{outbn}'):
            self.log.debug(f'already untarred {tarpath}- doing nothing')    
        else :                
            cmd = ['tar', '-C',os.path.dirname(tarpath), '-xvf', tarpath ]
            run_command(cmd)
            self.log.debug(f'untarred {tarpath}') 

        if not self.nocleanup:
            os.remove(tarpath)

        # return(stderr, stdout, returncode)

    def get_read_length(self, gzpath):
        if os.path.isfile(gzpath):
            z_head = subprocess.Popen( f'zcat {gzpath} | head',
                stdout=subprocess.PIPE, shell=True ).communicate()
            
            stdout = z_head[0].decode().split('\n')
            l = len(stdout[1])

            return(l)


def get_ids(fastqs = FASTQS.srcurl) :

    # assign project ids
    proj_ids = [os.path.dirname(f).replace('http://data.nemoarchive.org/biccn/grant/','') for f in fastqs  ]
    
    proj_ids = [re.sub('transcriptome/','',pid).replace('mouse/','').replace('raw/','') for pid in proj_ids ]
    proj_ids = [re.sub('/' ,'_',pid) for pid in proj_ids]

    # get the run id - separates lanes ideally.
    run_ids = [os.path.basename(f).replace('.fastq.tar','') for f in fastqs  ]
    run_ids = [re.sub('.fastq.tar','',r) for r in run_ids  ]
    run_ids = [re.sub('_L00[0-9]','',r) for r in run_ids  ]

    pdf = pd.DataFrame({'proj_id':proj_ids})
    pdf['run_id'] = run_ids
    pdf['tech'] = 'SSv4'
    pdf.loc[pdf.proj_id.str.upper().str.contains('10X_V2'), 'tech']  = '10x_v2'
    pdf.loc[pdf.proj_id.str.upper().str.contains('10X_V3'), 'tech']  = '10x_v3'

    # assign the smartseq 'run_ids' as the 'proj_id'. 
    pdf.loc[pdf.tech =='SSv4', 'run_id'] = pdf.loc[pdf.tech =='SSv4', 'proj_id']

    pdf['uid'] = [os.path.basename(f).replace('.fastq.tar','') for f in fastqs]
    pdf.to_csv('/home/johlee/git/scqc/resource/biccn_projects.tsv',sep="\t")
    return(pdf)

# get_ids()
    
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
