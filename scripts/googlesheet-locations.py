#!/usr/bin/env python
#
#  
#
#
#
import io
import pandas as pd
import requests

# google sheet URL 
url = "https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/export?gid=0&format=tsv"

r = requests.post(url) 
with io.StringIO(r.text) as imf:
    df = pd.read_csv(imf, sep='\t')
#print(df)
locations = list(df['Data location'][df['Data location'].notna()])
locations.sort()
for s in locations:
    print(f"{s}")