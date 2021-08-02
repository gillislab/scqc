#!/usr/bin/env python

# get all Supplement data - generally only needs to be done once
# getGenomeData() uses os.system(wget ...) to get data. FTP + gz links were annoying me.
# getSTAR() and getSRAtoolkit() don't do anything. Just comments on how to download them for now.


import gzip
import os

import requests

#### config -set up  #####
species = "mouse"
nCore = 30

# directories to store data
suppDirec = "SupplementData"
metaDirec = "MetaData" 
fastqDirec = "FASTQ"
staroutDirec ="STARout"
statsDirec = "QCStats"

# Don't change keys! Just the value to define the directory name
suppDirecSubs = {   
    "genomes" :"genomes",       # where to store genomes 
    "whitelists" :"whitelists", # where to store whitelists
    "starindex": "STARindices"  # where to store star indices
}


# Don't change keys! Just the value to define the directory name

# each subdirectory will contain a dataframe for each project in the given stage. 
# DFs are moved to the next stage upon completion
metaDirecSubs = {       # order matters!
    "query" : "Efetch_Queried" ,  # where to store efetch DFs
    "tp" : "Tech_Predicted",   # where to store LCP parsed DFS
    "fastq" :"FASTQ_Downloaded", # where to store FASTQ DFs
    "staraligned":"STAR_Aligned",     # where to store STAR_aligned DFs
    "datadone" :"Data_Saved"        # where to store Data_saved DFs
}

#################




def MakeDataDirecs(suppDirec="SupplementData", metaDirec ="MetaData", fastqDirec="FASTQ",staroutDirec="STARout",suppDirecSubs= ["genomes","whitelists"], metaDirecSubs=[    "Efetch_Queried" ,"Tech_Predicted","FASTQ_Downloaded","STAR_Aligned","Data_Saved"],statsDirec ="QCstats", species="mouse") : 
    '''
    Set up the directory. 

    '''

    for sub in suppDirecSubs :
        try: 
            os.makedirs("/".join([suppDirec, sub]) )
        except FileExistsError :
            print("Directory exists: '"+ suppDirec +"'. Doing nothing.")

    for sub in metaDirecSubs : 
        try: 
            dirname ="/".join([metaDirec,species,sub])
            os.makedirs(dirname  )
        except FileExistsError :
            print("Directory exists: '"+ dirname +"'. Doing nothing.")


    try: 
        os.makedirs("/".join([fastqDirec,species]) )
    except :
        pass

    try: 
        os.makedirs("/".join([staroutDirec,species]) )
    except FileExistsError :
        print("Directory exists: '"+ "/".join([fastqDirec,species]) +"'. Doing nothing.")

    try: 
        os.makedirs("/".join([statsDirec,species]) )
    except FileExistsError :
        print("Directory exists: '"+ "/".join([fastqDirec,species]) +"'. Doing nothing.")





    return 

def getWhiteLists(suppDirec ="SupplementData", whiteListSubDirec = "whitelists" ):  
    
    outDirec = "/".join([suppDirec, whiteListSubDirec])
    urls = {
        "whitelist_10xv1.txt" : "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt",
        "whitelist_10xv2.txt" : "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt",
        "whitelist_10xv3.txt" : "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
    }

    keys = list(urls.keys())

    for key in keys :

        url = urls[key]
        if not os.path.exists("/".join([outDirec, key])) :
            r = requests.get(url)
            
            if r.status_code == 200 :
                with open("/".join([outDirec, key]) ,"w" ) as f : 

                    if url.endswith(".gz") :
                        f.write(gzip.decompress(r.content).decode())
                    else :
                        f.write(r.text)
            else : 
                print( "Requesting "+ url +" failed with status_code: "+str(r.status_code))
            
    return

# currently uses os.system to get fasta/gtf files and decompress....
def getGenomeData(species, suppDirec = "SupplementData" , genomeSubDirec = "genomes") :
    outdir = "/".join([suppDirec, genomeSubDirec,species])
    try :
        os.makedirs( outdir)
    except FileExistsError :
        print("Directory exists: '"+ outdir +"'. Doing nothing.")



    spec2urls ={
        "human" :{
            "genome.fa":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz",
            "annotation.gtf":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz"
        },
        "mouse" :{
            "genome.fa":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.primary_assembly.genome.fa.gz",
            "annotation.gtf":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz"
        }
    }


    
    keys = list(spec2urls[species].keys())
    paths={}
    for key in keys :
        url = spec2urls[species][key]
        # does it already exist?
        if os.path.exists(outdir+"/"+key): 
            print( "File already exists. Doing nothing. '" + outdir+"/"+key+"'")
        else :
            os.system( "wget -nv -nc -O" + outdir+"/"+key+".gz" +" "+url  )
            os.system(" gunzip "+ outdir+"/"+key+".gz" )

        paths[key] = outdir+"/"+key    

    

    return (paths)

def generateGenomeIndices(species, suppDirec ="SupplementData", genomeSubDirec ="genomes",starindex = "STARindices",nCore = 30) :

    paths = getGenomeData(species ,suppDirec, genomeSubDirec)
    
    outdir = "/".join([suppDirec,genomeSubDirec,species  ,starindex])
    try:
        os.makedirs( outdir)
    except FileExistsError : 
        print("Directory exists: '"+ oudir+ "'. Doing nothing.")


    os.system(
        "STAR" + 
        " --runMode genomeGenerate " + 
        " --runThreadN " + str(nCore) +
        " --genomeDir " + outdir +
        " --genomeFastaFiles " + paths['genome.fa']  +
        " --sjdbGTFfile " +paths["annotation.gtf"]
    )
    

def getSRAtoolkit() : 
    # shell commands 

    # wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    # tar -vxzf sratoolkit.tar.gz
    # export PATH=$PATH:$PWD/sratoolkit*/bin
    pass

def getSTAR() : 
    # shell commads

    # wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
    # tar -xzf 2.7.8a.tar.gz
    # cd STAR-2.7.8a/source
    # make STAR
    # export PATH=$PATH:$PWD/STAR-2.7.8a/bin/Linux_x86_64
    pass


def main( 
        species = "mouse",
        nCore = 30,
        suppDirec = "SupplementData",
        metaDirec = "MetaData" ,
        fastqDirec = "FASTQ",
        staroutDirec ="STARout",
        statsDirec = "QCStats",
        suppDirecSubs = {   
            "genomes" :"genomes",       
            "whitelists" :"whitelists", 
            "starindex": "STARindices"  
        }, 
        metaDirecSubs = {       # order matters!
            "query" : "Efetch_Queried" ,  # where to store efetch DFs
            "tp" : "Tech_Predicted",   # where to store LCP parsed DFS
            "fastq" :"FASTQ_Downloaded", # where to store FASTQ DFs
            "staraligned":"STAR_Aligned",     # where to store STAR_aligned DFs
            "datadone" :"Data_Saved"        # where to store Data_saved DFs
        }
    ) :

    MakeDataDirecs(
        suppDirec, 
        metaDirec, 
        fastqDirec,
        staroutDirec,
        list(suppDirecSubs.values()), 
        list(metaDirecSubs.values()),
        statsDirec,
        species
    )


    getWhiteLists(
        suppDirec , 
        suppDirecSubs["whitelists"] 
    )

    generateGenomeIndices(
        species,
        suppDirec, 
        suppDirecSubs["genomes"],
        suppDirecSubs["startindex"],
        nCore 
    )