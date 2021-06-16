import gzip
import os
import logging
import shutil
import tempfile
import traceback
import urllib

from ftplib import FTP

import pandas as pd

def readlist(filepath):
    '''
    Assumes file is a list of strings, one per line. 
    '''

    if filepath is not None:
        logging.info(f'reading file: {filepath}')
        flist = []
        try:
            with open(filepath, 'r') as f:
                flist = [line.strip() for line in f]
            logging.debug(f'got list with {len(flist)} items.')
            return flist
        except:
            return []
    else:
        logging.info('no file. return [].')
        return []


def writelist(filepath, dlist):
    logging.info(f"writing list length={len(dlist)} to file='{filepath}'")
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        with os.fdopen(tfd, 'w') as f:
            nlines = 0
            for item in dlist:
                f.write(f"{item}\n")
                nlines += 1
        os.rename(tfname, filepath)
        logging.info(f"wrote {nlines} to {filepath}")
    except Exception as ex:
        logging.error(traceback.format_exc(None))

    finally:
        pass

def merge_write_df(newdf, filepath):
    """
    Reads existing, merges new, drops duplicates, writes to temp, renames temp. 
    """
    log = logging.getLogger('utils')
    log.debug(f'inbound new df: {newdf}')
    if os.path.isfile(filepath):
        df = pd.read_csv(filepath, sep='\t', index_col=0, comment="#")
        log.debug(f'read df: {df}')
        df = df.append(newdf, ignore_index=True )
        log.debug(f'appended df: {df}')
    else:
        df = newdf
    df.drop_duplicates(inplace=True)

    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        df.to_csv(tfname, sep='\t')
        
        
        os.rename(tfname, filepath)
        logging.info(f"wrote df to {filepath}")
   
    except Exception as ex:
        logging.error(traceback.format_exc(None))



def listdiff(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 - s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"diff has length {len(dl)}")
    return dl


def listmerge(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 | s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"merged has length {len(dl)}")
    return dl


def download_ftpurl(srcurl, destpath, finalname=None, overwrite=True, decompress=True):
    """
    Downloads via FTP from ftp src url to local destpath, 
    If finalname is specified, renames final output. 
    overwrite: won't re-download if filename already exists. 
    decompress: if filename ends with .gz , will gunzip  
    """
    log = logging.getLogger('star')
    (scheme, host, fullpath, p, q ,f) = urllib.parse.urlparse(srcurl)
    filename = os.path.basename(fullpath)
    dirname = os.path.dirname(fullpath)    
    log.info(f"Downloading file {filename} at path {dirname}/ on host {host} via FTP.")
    ftp = FTP(host)
    ftp.login('anonymous','hover@cshl.edu')
    ftp.cwd(dirname)
    log.debug(f'opening file {destpath}/{filename}. transferring...')
    with open(f'{destpath}/{filename}', 'wb') as fp:
        ftp.retrbinary(f'RETR {filename}', fp.write)
    log.debug(f"done retrieving {destpath}/{filename}")
    ftp.quit()
    
    if decompress and filename.endswith('.gz'):
        log.debug(f'decompressing gzip file {destpath}/{filename}')
        gzip_decompress(f'{destpath}/{filename}')
        os.remove(f'{destpath}/{filename}')
        filename = filename[:-3]
    
    if finalname is not None:
        src = "/".join([ destpath , filename])
        dest = "/".join([ destpath , finalname])
        log.info(f'renaming {src} -> {dest}')
        os.rename(src, dest )

    
def gzip_decompress(filename):
    """
    default for copyfileobj is 16384
    https://blogs.blumetech.com/blumetechs-tech-blog/2011/05/faster-python-file-copy.html
    
    """
    log = logging.getLogger('utils')
    if filename.endswith('.gz'):
        targetname = filename[:-3]
        bufferlength = 10 * 1024 * 1024 # 10 MB
        with gzip.open(filename, 'rb') as f_in:
            with open(targetname, 'wb') as f_out:
                shutil.copyfileobj(f_in,f_out, length=bufferlength)
    else:
        log.warn(f'tried to gunzip file without .gz extension {filename}. doing nothing.')
