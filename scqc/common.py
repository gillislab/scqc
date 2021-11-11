
import re



# Translate between Python and SRAToolkit log levels for wrapped commands.
#  fatal|sys|int|err|warn|info|debug
LOGLEVELS = {
    10: 'debug',
    20: 'info',
    30: 'warn',
    40: 'err',
    50: 'fatal',
}


PROJ_COLUMNS = ['proj_id', 'ext_ids', 'title', 'abstract', 'submission_id', 'data_source']

SAMP_COLUMNS = ['samp_id', 'ext_ids',  'taxon',
                'sciname', 'title', 'attributes', 'proj_id', 'submission_id', 'data_source']

EXP_COLUMNS = ['exp_id', 'ext_ids',  'strategy',
               'source', 'lcp', 'samp_id', 'proj_id', 'submission_id', 'data_source']

RUN_COLUMNS = ['run_id', 'ext_ids', 'tot_spots', 'tot_bases', 'run_size', 'publish_date',
               'taxon', 'organism', 'nreads',  'basecounts', 'file_url','file_size','exp_id', 'samp_id', 'proj_id', 
               'submission_id', 'data_source' ]

IMPUTE_COLUMNS = ['run_id' ,'tech_version','read1','read2','exp_id','samp_id','proj_id', 'taxon', 'batch', 'data_source']


TECH_RES = {
    '10x'   : re.compile("10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'|Kit v1|PN-120233|10X V1", re.IGNORECASE),
    #'10xv1' : re.compile("", re.IGNORECASE),
    #'10xv2' : re.compile("v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2", re.IGNORECASE),
    #'10xv3' : re.compile("v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078", re.IGNORECASE),
    'smartseq' : re.compile("Smart-Seq|SmartSeq|Picelli|SMART Seq", re.IGNORECASE),
    'smarter' : re.compile("SMARTer", re.IGNORECASE),
    'dropseq' : re.compile("Cell 161, 1202-1214|Macosko|dropseq|drop-seq", re.IGNORECASE),
    'celseq'  : re.compile("CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq", re.IGNORECASE),
    'sortseq' : re.compile("Sort-seq|Sortseq|Sort seq", re.IGNORECASE),
    'seqwell' : re.compile("Seq-Well|seqwell", re.IGNORECASE),
    'biorad'  : re.compile("Bio-Rad|ddSeq", re.IGNORECASE),
    'indrops' : re.compile("inDrop|Klein|Zilionis", re.IGNORECASE),
    'marsseq2': re.compile("MARS-seq|MARSseq|Jaitin et al|jaitin", re.IGNORECASE),
    'tang'    : re.compile("Tang", re.IGNORECASE),
    # 'TruSeq':re.compile("TruSeq", re.IGNORECASE),
    'splitseq': re.compile("SPLiT-seq", re.IGNORECASE),
    'microwellseq': re.compile("Microwell-seq", re.IGNORECASE)
}

NEMO_URL_TECH_MAP = { 'SSv4' : 'smartseq',
                '10x_v1': '10xv1',
                '10x_v2': '10xv2',                
                '10x_v3': '10xv3',
                }


#deprecated
keywords = {
    "is10x": "10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'",
    "v3": "v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078",
    "v2": "v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2",
    "v1": "Kit v1|PN-120233|10X V1",
    "ss": "Smart-Seq|SmartSeq|Picelli|SMART Seq",
    "smarter": "SMARTer",
    "dropseq": "Cell 161, 1202-1214|Macosko|dropseq|drop-seq",
    "celseq": "CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq",
    "sortseq": "Sort-seq|Sortseq|Sort seq",
    "seqwell": "Seq-Well|seqwell",
    "biorad": "Bio-Rad|ddSeq",
    "indrops": "inDrop|Klein|Zilionis",
    "marsseq2": "MARS-seq|MARSseq|Jaitin et al|jaitin",
    "tang": "Tang",
    # "TruSeq":"TruSeq",
    "splitseq": "SPLiT-seq",
    "microwellseq": "Microwell-seq"
}




def get_project_backend(df, proj_id):
    list(df[df.proj_id == proj_id]['data_source'])[0]