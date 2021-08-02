#!/usr/bin/env python

import sys

print(f'python ver. {sys.version} ')

from urllib import parse

db='sra'
#species = ['"homo sapiens"[Organism]','"mus musculus"[Organism]']
species = ['"mus musculus"[Organism]']
strategy = ['"rna seq"[Strategy]']
textword = ['"single cell"[Text Word]','"brain"[Text Word]']

u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra" 
spc = ' OR '.join(species)
print(f'spc = {spc}')
strat = ' OR '.join(strategy)
print(f'strat = {strat}')
txtwrd = ' AND '.join(textword)
print(f'txtwrd = {txtwrd}')
tis = '"tissue"[Attribute Name]'

params = f"({spc}) AND ({strat}) AND ({txtwrd}) "
print(f'uncoded params = {params}')

encoded = parse.quote_plus(params)
furl = f"{u}&term={encoded}&retmode=xml"

print(furl)
