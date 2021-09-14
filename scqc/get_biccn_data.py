import os 
import glob
import sys
import subprocess
import argparse
import requests
import pandas as pd
import logging
import datetime as dt
import time

# search request
# https://portal.nemoarchive.org/search/s?facetTab=cases&filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.modality%22,%22value%22:%5B%22Transcriptomics%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.organism%22,%22value%22:%5B%22Mouse%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.subspecimen_type%22,%22value%22:%5B%22Nuclei%22,%22Cells%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22file.format%22,%22value%22:%5B%22FASTQ%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.technique%22,%22value%22:%5B%22Patch-seq;SMART-seq%20v4%22,%22SMART-seq%20v4%22,%2210x%20Chromium%203%27%20v3%20sequencing%22,%2210x%20Chromium%203%27%20v2%20sequencing%22%5D%7D%7D%5D%7D

os.chdir('/data/biccn/nemo/')

class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """


def download_wget(srcurl, destpath, finalname=None, overwrite=True, decompress=True, rate='1M'):
    """
    
    
    
    GNU Wget 1.20.1, a non-interactive network retriever.
    Usage: wget [OPTION]... [URL]...
    
    Startup:
      -V,  --version                   display the version of Wget and exit
      -h,  --help                      print this help
      -v,  --verbose                   be verbose (this is the default)
      -nv, --no-verbose                turn off verboseness, without being quiet
           --report-speed=TYPE         output bandwidth as TYPE.  TYPE can be bits
      -t,  --tries=NUMBER              set number of retries to NUMBER (0 unlimits)
           --retry-connrefused         retry even if connection is refused
           --retry-on-http-error=ERRORS    comma-separated list of HTTP errors to retry
      -O,  --output-document=FILE      write documents to FILE
      -nc, --no-clobber                skip downloads that would download to
                                         existing files (overwriting them)
    
      -c,  --continue                  resume getting a partially-downloaded file
           --progress=TYPE             select progress gauge type
           --show-progress             display the progress bar in any verbosity mode
      -N,  --timestamping              don't re-retrieve files unless newer than
                                         local
           --no-if-modified-since      don't use conditional if-modified-since get
                                         requests in timestamping mode
           --no-use-server-timestamps  don't set the local file's timestamp by
                                         the one on the server
       -T,  --timeout=SECONDS           set all timeout values to SECONDS
           --dns-timeout=SECS          set the DNS lookup timeout to SECS
           --connect-timeout=SECS      set the connect timeout to SECS
           --read-timeout=SECS         set the read timeout to SECS
      -w,  --wait=SECONDS              wait SECONDS between retrievals
           --waitretry=SECONDS         wait 1..SECONDS between retries of a retrieval
           --random-wait               wait from 0.5*WAIT...1.5*WAIT secs between retrievals
    
           --limit-rate=RATE           limit download rate e.g. 1M  1 MB/s      
    """
    logging.debug(f'wget file {srcurl}')
    cmd = ['wget',
        #    '--no-verbose',
           '--no-use-server-timestamps',
           '--limit-rate', rate,
           '--continue', 
           '-O', f'{destpath}',
           f'{srcurl}']
    cmdstr = " ".join(cmd)
    logging.debug(f"wget command: {cmdstr} running...")
    
    start = dt.datetime.now()
    cp = subprocess.run(cmd, 
                        universal_newlines=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == '0':
        logging.debug(f"got stderr: {cp.stderr}")
        logging.debug(f"got stdout: {cp.stdout}")
        if len(cp.stderr) > 10:
            pass
            # dlbytes = parse_wget_output_bytes(cp.stderr)
            # logging.info(f'downloaded {dlbytes} bytes {destpath} successfully, in {elapsed.seconds} seconds. ')
        else:
            logging.info(f'file already downloaded.')
    else:
        logging.error(f'non-zero return code for src {srcurl}')
    return cp.returncode

def parse_wget_output_bytes(outstr):
    """
    E.g. 2021-07-20 14:33:09 URL:https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5529542/SRR5529542 [17019750/17019750] -> "SRR5529542.sra" [1]
    """
    logging.debug(f'handling stderr string {outstr}')    
    fields = outstr.split()
    bstr = fields[3][1:-1]
    dlbytes = int(bstr.split('/')[0])
    return dlbytes

def run_command(cmd):
    """
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.info(f"command: {cmdstr} running...")
    start = dt.datetime.now()
    cp = subprocess.run(cmd, 
                    text=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.debug(f"got stderr: {cp.stderr}")
        print(cp.stderr)
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")
    
    if str(cp.returncode) == '0':
        logging.info(f'successfully ran {cmdstr}')
    else:
        print(f'non-zero return code for cmd {cmdstr}')
        
        # raise NonZeroReturnException()

    return(cp.stderr, cp.stdout ,cp.returncode)

def get_files_from_html(htmlpath,pattern = '/'):
    # note pattern = "/" looks for subdirectories only.

    with open(htmlpath) as f :
        htmltext = f.read().split('<a href="')
        subdirs = [ txt.split('</a>')[0] for txt in htmltext]
        subdirs = [ subdir.split('">')[1] for subdir in subdirs ][1:]
        subdirs = [ subdir  for subdir in subdirs if pattern in subdir]
    return(subdirs)
                


def get_zeng_data(base_url = 'http://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome'):

    subsamples = ['scell','sncell']
    techs = ['10x_v2','10x_v3','SSv4','SSv4_viral']
    species = ['mouse']
    allfastqs = []
    # note : sncell/SSv4_viral does not exist. failure expected
    for spec in species :
        for subsamp in subsamples:
            for tech in techs: 
                srcurl = f'{base_url}/{subsamp}/{tech}/{spec}/raw'
                
                # first, get the index.html file
                cmd = ['wget',
                    '--recursive',
                    '--level','1',
                    f'{srcurl}/'
                ]   
                stderr, stdout, returncode =run_command(cmd)
                
                if str(returncode) == '0' :
                    # next, parse the index.html to get the tissue names.
                    htmlpath =(f'{srcurl}/index.html').replace('http://','')
                    regions = get_files_from_html(htmlpath,pattern ="/")

                    for region in regions :
                        cmd = ['wget',
                            '--recursive',
                            '--level','1',
                            f'{srcurl}/{region}'
                        ]
                        stderr, stdout ,returncode =run_command(cmd)
                        if str(returncode) =='0':
                            htmlpath =(f'{srcurl}/{region}index.html').replace('http://','')
                            fastqs = get_files_from_html(htmlpath,pattern = 'fastq')
                            fastqs = [ f'{srcurl}/{region}{fastq}' for fastq in fastqs  ] 
                            allfastqs = allfastqs + fastqs
                else: 
                    print(f'something wrong with {srcurl}/{region}')

    
    f = open('/home/johlee/zeng_fastqs.txt','w')
    [f.write( f'{fq}\n' ) for fq in allfastqs ]
    f.close()
    return(allfastqs)



                        
    # scell/10x_v2/mouse/raw/*
    # scell/10x_v3/mouse/raw/MOp/
    # scell/SSv4/mouse/raw/*
    # scell/SSv4_viral/mouse/raw/*

    # sncell/10x_v2/mouse/raw/MOp
    # sncell/10x_v3/mouse/raw/MOp
    # sncell/SSv4/mouse/raw/MOp/

def get_huang_data(base_url = 'http://data.nemoarchive.org/biccn/grant/u19_huang'):
    pis = ['arlotta','dulac','macosko','macosko_regev'] 
    modals = ['transcriptome']
    subsamples = ['scell','sncell']
    techs = ['10x_v2','10x_v3']

    allfastqs = []
    for pi in pis:
        for modal in modals:
            for subsamp in subsamples:
                for tech in techs:
                    if pi =='dulac':
                        srcurl = f'{base_url}/{pi}/{modal}/{subsamp}/{tech}/mouse/pag/raw/'
                    elif pi =='arlotta' :               
                        srcurl = f'{base_url}/{pi}/{modal}/{subsamp}/{tech}/mouse/raw/'
     
                    if pi == 'dulac' or pi =='arlotta':
                        cmd = ['wget',
                            '--recursive',
                            '--level','1',
                            f'{srcurl}'
                        ]   
                        stderr, stdout, returncode =run_command(cmd)
                        # index html exists for dulac
                        if str(returncode) == '0':
                            htmlpath =(f'{srcurl}/index.html').replace('http://','')
                            fastqs = get_files_from_html(htmlpath,pattern = 'fastq')

                            fastqs = [ f'{srcurl}{fq}'  for fq in fastqs] 
                            allfastqs += fastqs

                    else : # macosko et al
                        tech2 = tech[:2] + tech[2].upper() + tech[3:]
                        # first look for the subdirecs 
                        srcurl = f'{base_url}/{pi}/{modal}/{subsamp}/{tech2}/mouse/raw/'
                        cmd = ['wget',
                            '--recursive',
                            '--level','1',
                            f'{srcurl}'
                        ]   
                        stderr, stdout, returncode =run_command(cmd)
                        if str(returncode) == '0':
                            htmlpath =(f'{srcurl}index.html').replace('http://','')
                            subdirs = get_files_from_html(htmlpath,pattern = '/')

                            for sd in subdirs :
                                cmd = ['wget',
                                    '--recursive',
                                    '--level','1',
                                    f'{srcurl}{sd}'
                                ]   
                                stderr, stdout, returncode =run_command(cmd)
                                htmlpath =(f'{srcurl}{sd}index.html').replace('http://','')
                                fastqs = get_files_from_html(htmlpath,pattern = 'fastq')

                                fastqs = [ f'{srcurl}{fq}'  for fq in fastqs] 
                                allfastqs += fastqs

            
    f = open('/home/johlee/huang_fastqs2.txt','w')
    [f.write( f'{fq}\n' ) for fq in allfastqs ]
    f.close()

    return(allfastqs)

def get_nowakowski_data(base_url ='http://data.nemoarchive.org/biccn/grant/rf1_nowakowski/nowakowski/transcriptome/scell/10x_v3/mouse') :
    regions = ['PNdev', 'evobc']
    
    allfastqs=[]
    for region in regions :
        srcurl = f'{base_url}/{region}/raw/'
        # get index htmls 
        cmd = ['wget',
            '--recursive',
            '--level','1',
            f'{srcurl}'
        ]   
        
        stderr, stdout, returncode =run_command(cmd)
        # index html exists for dulac
        if str(returncode) == '0':
            htmlpath =(f'{srcurl}/index.html').replace('http://','')
            fastqs = get_files_from_html(htmlpath,pattern = 'fastq')

            fastqs = [ f'{srcurl}{fq}'  for fq in fastqs] 
            allfastqs += fastqs

    f = open('/home/johlee/nowakowski_fastqs.txt','w')
    [f.write( f'{fq}\n' ) for fq in allfastqs ]
    f.close()

    return(allfastqs)


df = pd.read_csv('/home/johlee/biccn_fastq_paths.txt',header=None)

# for f in df[0][::-1]:
#     print(f)
#     rc = download_wget(f,destpath=f'/data/biccn/nemo/{os.path.basename(f)}')
#     time.sleep(1)
