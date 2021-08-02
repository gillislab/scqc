# pick a project
import pandas as pd
import os
species = "mouse"   # mouse or human for now           

getAll = False          # Boolean  
getRecent = True        # Boolean - ignored if getAll == True
lastNdays = 5           # only used if getRecent == True
nCore = 30              # max cores allows for this project accession
suppDirec = "/home/johlee/SCQC/SupplementData"    # supplemental data directory. 
metaDir = "/home/johlee/SCQC/MetaData"            # MetaData Directory - should contain <run>_MetaData_TP.tsv 
fastqStore="/home/johlee/SCQC/FASTQ"              # Temporarily store fastq files here
STARout = "/home/johlee/SCQC/STARout"             # Temporily store star outputs 


project_acc = 'SRP258272' 



fastqstoredir = "/".join([fastqStore, project_acc])
def tax2spec(taxid):
    d= {    "10090" :"mouse",
            "9606" : "human"}
    return(d[str(taxid)])


def getSTARParameters(version = "10xv3",suppDirec ="SupplementData") :
    d = {
        "10xv1": 
            {        
                "solo_type":"CB_UMI_Simple",
                "whiteListPath" : suppDirec+"/whitelists/whitelist_10xv1.txt",
                "CB_length":"14",
                "UMI_start":"15",
                "UMI_length":"10" 
            },
        "10xv2":
            {        
                "solo_type":"CB_UMI_Simple",
                "whiteListPath" : suppDirec+"/whitelists/whitelist_10xv2.txt",
                "CB_length":"16",
                "UMI_start":"17",
                "UMI_length":"10" 
            },
        "10xv3" : 
            {        
                "solo_type":"CB_UMI_Simple",
                "whiteListPath" : suppDirec+"/whitelists/whitelist_10xv3.txt",
                "CB_length":"16",
                "UMI_start":"17",
                "UMI_length":"12" 
            },
        "SMART-Seq" :{
                "solo_type":"SmartSeq",
                "soloUMIdedup":"Exact",
                "soloStrand" :"Unstranded"
        }
    }
    return(d[version])



def align_10x(project_acc  ,df, metaDir ,suppDirec, STARout , nCore=30):
    filename = "/".join([metaDir, project_acc]) + "_10x_MetaData.tsv"
    df = pd.read_csv(filename, sep="\t")

    for i in  range(len(df.Run)) : 

        spec = tax2spec(df.Taxon_ID.iloc[i])
        readsIn =df.readsIn.iloc[i] 
        method = df.Method.iloc[i] 
        run = df.Run.iloc[i] 
        param = get10xParameters(method , suppDirec)
        genomeDir = "/".join([suppDirec , "genomes",spec,"STARindices"])

        os.system( (
            "STAR" + 
            " --runMode alignReads" +
            " --runThreadN "+ str(nCore)+
            " --genomeDir "+ genomeDir+
            " --outFileNamePrefix "+ STARout+ "/" + run  +
            " --soloType " +param['solo_type']+
            " --soloCBwhitelist "+ param['whiteListPath']+
            " --soloCBlen "+ param["CB_length"]+
            " --soloUMIlen " +param["UMI_length"]+
            " --soloUMIstart "+ param["UMI_start"]+
            " --soloFeatures Gene "+
            # " --readFilesPrefix " + fastqstoredir +"/" +run +
            " --readFilesIn "+  readsIn
            # " --readFilesManifest ", paste0(path2fastq,"/",manifest)
            )
        )


def align_SS(project_acc  ,df, metaDir ,suppDirec, STARout , nCore=30):
    filename = "/".join([metaDir, project_acc]) + "_SS_MetaData.tsv"
    df = pd.read_csv(filename, sep="\t")

    # for i in  range(len(df.Run)) : 

    spec = tax2spec(df.Taxon_ID.unique)
    readsIn =df.readsIn.unique
    method = df.Method.unique
    run = df.Run.unique
    param = getSTARParameters(method , suppDirec)
    genomeDir = "/".join([suppDirec , "genomes",spec,"STARindices"])
    manifestPath = ""
    
    os.system( (
        "STAR"  
        + " --runMode alignReads" 
        + " --runThreadN "+ str(nCore)
        + " --genomeDir "+ genomeDir
        + " --outFileNamePrefix "+ STARout+ "/" + run  
        + " --soloType " +param['solo_type']
        + " --soloUMIdedup " + param["soloUMIdedup" ] 
        + " --soloStrand " + param["soloStrand"] 
        + " --soloFeatures Gene "
        + " --readFilesManifest ", manifestPath
        )
    )



def main(): 
    
    pass