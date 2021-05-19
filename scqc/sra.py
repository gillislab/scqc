#!/usr/bin/env python
#
#
# http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=
#
#

import io
import logging
import pandas
import requests


def get_run_metadata(sraproject):
    '''
    
    
    '''
    url="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi"

    headers = {
        "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding":"gzip,deflate,sdch",
        "Accept-Language":"en-US,en;q=0.8",
        "Cache-Control":"no-cache",
        "Connection":"keep-alive",
        "DNT":"1",
        "Host":"trace.ncbi.nlm.nih.gov",
        "Origin":"http://trace.ncbi.nlm.nih.gov",
        "Pragma":"no-cache",
        "Referer":"http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?db=sra",
        "User-Agent":"Mozilla/5.0 (iPhone; CPU iPhone OS 6_0 like Mac OS X) AppleWebKit/536.26 (KHTML, like Gecko) Version/6.0 Mobile/10A5376e Safari/8536.25"}
    
    payload = {
        "db":"sra",
        "rettype":"runinfo",
        "save":"efetch",
        "term": sraproject }
    
    r = requests.put(url, data=payload, headers=headers, stream=True) 
    with io.BytesIO(r.content) as imf:
        df = pandas.read_csv(imf)
    return df


def get_run_files(srarun):
    '''
    
    '''
    pass
    #response = requests.get(url, stream=True)
    #handle = open(target_path, "wb")
    #for chunk in response.iter_content(chunk_size=512):
    #    if chunk:  # filter out keep-alive new chunks
    #        handle.write(chunk)


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.DEBUG)

    SRAPROJ="SRP185852"
    df = get_run_metadata(SRAPROJ)
    logging.debug(f"Got list of {len(df)} runs")
    runs = list(df['Run'])
    logging.info(f"Runlist: {runs}")


