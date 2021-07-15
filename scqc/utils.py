import gzip
import os
import logging
import shutil
import subprocess
import tempfile
import traceback
import urllib
import numpy as np
from scipy import sparse, stats
from ftplib import FTP
from configparser import ConfigParser
import io
import bottleneck
import pandas as pd


def readlist(filepath):
    '''
    Assumes file is a list of strings, one per line. 
    Ignores lines beginning with a has '#'
    Ignores characters in a line afeter a '#'
    '''

    if filepath is not None:
        logging.info(f'reading file: {filepath}')
        flist = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    idx = line.find('#')
                    if idx == -1:
                        flist.append(line.strip())
                    elif idx > 0:
                        flist.append(line[:idx].strip())
                        
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


def load_df(filepath):
    """
    Convenience method to load DF
    """
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    return df



def merge_write_df(newdf, filepath):
    """
    Reads existing, merges new, drops duplicates, writes to temp, renames temp. 
    """
    log = logging.getLogger('utils')
    log.debug(f'inbound new df: {newdf}')
    if os.path.isfile(filepath):
        df = pd.read_csv(filepath, sep='\t', index_col=0) #, comment="#"
        log.debug(f'read df: {df}')
        df = df.append(newdf, ignore_index=True)
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
           '--limit-rate', rate,
           '--continue', 
           '-O', f'{destpath}',
           f'{srcurl}']
    cmdstr = " ".join(cmd)
    logging.debug(f"wget command: {cmdstr} running...")
    cp = subprocess.run(cmd)
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == '0':
        logging.info(f'downloaded {destpath} successfully.')
    else:
        self.log.error(f'non-zero return code for src {srcurl}')
    return cp.returncode
    


def download_ftpurl(srcurl, destpath, finalname=None, overwrite=True, decompress=True):
    """
    destpath is directory
    
    Downloads via FTP from ftp src url to local destpath, 
    If finalname is specified, renames final output. 
    overwrite: won't re-download if filename already exists. 
    decompress: if filename ends with .gz , will gunzip  
    """
    log = logging.getLogger('star')
    # source FTP info
    (scheme, host, fullpath, p, q, f) = urllib.parse.urlparse(srcurl)
    filename = os.path.basename(fullpath)
    dirname = os.path.dirname(fullpath)
    
    # local files?
    localfile = f'{destpath}/{filename}'
    localfinal = f'{destpath}/{finalname}'
    destexists = os.path.exists(localfile) or os.path.exists(localfinal)
    log.debug(f'checking if {localfile} or {localfinal} exist -> {destexists}')
    
    if destexists and not overwrite:
        log.info(f"Destination files already exist and overwrite=false. Skipping.")
    else:
        log.info(
            f"Downloading file {filename} at path {dirname}/ on host {host} via FTP.")
        ftp = FTP(host)
        ftp.login('anonymous', 'hover@cshl.edu')
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
            src = "/".join([destpath, filename])
            dest = "/".join([destpath, finalname])
            log.info(f'renaming {src} -> {dest}')
            os.rename(src, dest)


def gzip_decompress(filename):
    """
    default for copyfileobj is 16384
    https://blogs.blumetech.com/blumetechs-tech-blog/2011/05/faster-python-file-copy.html

    """
    log = logging.getLogger('utils')
    if filename.endswith('.gz'):
        targetname = filename[:-3]
        bufferlength = 10 * 1024 * 1024  # 10 MB
        with gzip.open(filename, 'rb') as f_in:
            with open(targetname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out, length=bufferlength)
    else:
        log.warn(
            f'tried to gunzip file without .gz extension {filename}. doing nothing.')



def gini_coefficient_fast(X):
    """ 
        expects a CSR sparse matrix
        Sorting is O(n log n ) (here n is at most number of genes)
        loops over cells (m) instead of gene pairs. 
        Overall, at most O( m n log n)  but realistically, 
        density of 10% -> `effective n` is 0.1 * n
        
        only looks at nonzero elements
        Cells with no expression get a gini score of 0       
    """    
    # x = np.asarray(x)
    g = np.zeros(X.shape[0])    # ncells
    n = X.shape[1]          # ngenes
    for i in range(X.shape[0]): # loops for all cells
        # take the nonzero elements of the ith cell
        x = X[i,:]  
        x= x[:, x.indices].A.flatten()

        sorted_x = np.sort(x)   
        cumx = np.cumsum(sorted_x, dtype=float)

        if len(cumx) == 0 : # cell with zero expression - perfect equilibrium
            g[i] = 0
        else :
            g[i] =(n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return g




def sparse_pairwise_corr(A, B=None):
    """
    Compute pairwise correlation for sparse matrices. 
    Currently only implements pearson correlation.
    If A is N x P
       B is M x P
    Result is N+M x N+M symmetric matrix
    with off diagonal blocks as correlations between 
        elements in A with elements in B
    and main diagonal blocks as correlations between
        elements in A (or B) with elements in A (or B)
    """

    if B is not None:
        A = sparse.vstack((A, B), format='csr')

    A = A.astype(np.float64)
    n = A.shape[1]

    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)

    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))

    return coeffs



def taxon_to_spec(taxid= '10090'):
    d = {   '10090': "mouse",
            '9606':"human"}
    return(d[taxid])


def spec_to_taxon(spec="mouse"):
    d = {   "mouse":"10090",
            "human":"9606"}         
    return(d[spec])


# EGAD functions compliments of Ben
def rank(data, nan_val=.5):
    """Rank normalize data
    
    Rank standardize inplace 
    Ignores Nans and replace with .5
    
    Does not return 
    Arguments:
        data {np.array} -- Array of data
    
    """
    finite = np.isfinite(data)
    ranks = bottleneck.rankdata(data[finite]).astype(data.dtype)

    ranks -= 1
    top = np.max(ranks)
    ranks /= top
    data[...] = nan_val
    data[np.where(finite)] = ranks

    return(data)


def run_egad(go, nw, **kwargs):
    """EGAD running function
    
    Wrapper to lower level functions for EGAD

    EGAD measures modularity of gene lists in co-expression networks. 

    This was translated from the MATLAB version, which does tiled Cross Validation
    
    The useful kwargs are:
    int - nFold : Number of CV folds to do, default is 3, 
    int - {min,max}_count : limits for number of terms in each gene list, these are exclusive values


    Arguments:
        go {pd.DataFrame} -- dataframe of genes x terms of values [0,1], where 1 is included in gene lists
        nw {pd.DataFrame} -- dataframe of co-expression network, genes x genes
        **kwargs 
    
    Returns:
        pd.DataFrame -- dataframe of terms x metrics where the metrics are 
        ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    """
    assert nw.shape[0] == nw.shape[1] , 'Network is not square'
    assert np.all(nw.index == nw.columns) , 'Network index and columns are not in the same order'
    nw_mask = nw.isna().sum(axis=1) != nw.shape[1]
    nw = nw.loc[nw_mask, nw_mask].astype(float)
    np.fill_diagonal(nw.values, 1)
    return _runNV(go, nw, **kwargs)


def _runNV(go, nw, nFold=3, min_count=20, max_count=1000):

    #Make sure genes are same in go and nw
    genes_intersect = go.index.intersection(nw.index)

    go = go.loc[genes_intersect, :]
    nw = nw.loc[genes_intersect, genes_intersect]

    #Make sure there aren't duplicates
    duplicates = nw.index.duplicated(keep='first')
    nw = nw.loc[~duplicates, ~duplicates]

    go = go.loc[:, (go.sum(axis=0) > min_count) & (go.sum(axis=0) < max_count)]
    go = go.loc[~go.index.duplicated(keep='first'), :]

    roc = _new_egad(go.values, nw.values, nFold)

    col_names = ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    #Put output in dataframe
    return pd.DataFrame(dict(zip(col_names, roc)), index=go.columns)


def _new_egad(go, nw, nFold):

    #Build Cross validated Positive
    x, y = np.where(go)
    cvgo = {}
    for i in np.arange(nFold):
        a = x[i::nFold]
        b = y[i::nFold]
        dat = np.ones_like(a)
        mask = sparse.coo_matrix((dat, (a, b)), shape=go.shape)
        cvgo[i] = go - mask.toarray()
        
    CVgo = np.concatenate(list(cvgo.values()), axis=1)

    sumin = np.matmul(nw.T, CVgo)

    degree = np.sum(nw, axis=0)

    predicts = sumin / degree[:, None]

    np.place(predicts, CVgo > 0, np.nan)

    #Calculate ranks of positives
    rank_abs = lambda x: stats.rankdata(np.abs(x))
    predicts2 = np.apply_along_axis(rank_abs, 0, predicts)

    #Masking Nans that were ranked (how tiedrank works in matlab)
    predicts2[np.isnan(predicts)] = np.nan

    filtering = np.tile(go, nFold)

    #negatives :filtering == 0
    #Sets Ranks of negatives to 0
    np.place(predicts2, filtering == 0, 0)

    #Sum of ranks for each prediction
    p = np.nansum(predicts2, axis=0)

    #Number of predictions
    #Number of 1's masked for each GO term for each CV
    n_p = np.sum(filtering, axis=0) - np.sum(CVgo, axis=0)

    #Number of negatives
    #Number of GO terms - number of postiive
    n_n = filtering.shape[0] - np.sum(filtering, axis=0)
    roc = (p / n_p - (n_p + 1) / 2) / n_n
    U = roc * n_p * n_n
    Z = (np.abs(U - (n_p * n_n / 2))) / np.sqrt(n_p * n_n *
                                                (n_p + n_n + 1) / 12)
    roc = roc.reshape(nFold, go.shape[1])
    Z = Z.reshape(nFold, go.shape[1])
    #Stouffer Z method
    Z = np.nansum(Z, axis=0) / np.sqrt(nFold)
    #Calc ROC of Neighbor Voting
    roc = np.nanmean(roc, axis=0)
    P = stats.norm.sf(Z)

    #Average degree for nodes in each go term
    avg_degree = degree.dot(go) / np.sum(go, axis=0)

    #Calc null auc for degree
    ranks = np.tile(stats.rankdata(degree), (go.shape[1], 1)).T

    np.place(ranks, go == 0, 0)

    n_p = np.nansum(go, axis=0)
    nn = go.shape[0] - n_p
    p = np.nansum(ranks, axis=0)

    roc_null = (p / n_p - ((n_p + 1) / 2)) / nn

    return roc, avg_degree, roc_null, P


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
    return cp


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()