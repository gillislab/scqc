# get FASTQ files
# import pickle
import glob
import os

import numpy as np
import pandas as pd

# part 3: Get the fastq files for the specified project ID, verify technology
# downloading the fastq files is the significant bottleneck of the whole process... 

### Config ####


species = "mouse"   # mouse or human for now           

getAll = False          # Boolean 
getRecent = True        # Boolean - ignored if getAll == True
lastNdays = 5           # only used if getRecent == True
nCore = 30              # max cores allows for this project accession
suppDirec = "/home/johlee/SCQC/SupplementData"    # supplemental data directory. 
metaDirec = "/home/johlee/SCQC/MetaData"    # MetaData Directory - should contain <run>_MetaData_TP.tsv 
fastqStore="/home/johlee/SCQC/FASTQ"              # Temporarily store fastq files here
STARout = "/home/johlee/SCQC/STARout"             # Temporily store star outputs 

project_acc = "SRP258272"
# project_acc = 'SRP300819'       # 10x example
# project_acc = 'SRP288492'       # SS example
# project_acc = 'ERP122986'       # SS example that breaks
        # - failed to resolve accession 'ERR4354523' - no data ( 404 )
        # 2021-05-19T19:42:11 fasterq-dump.2.11.0 err: invalid accession 'ERR4354523'
        # fasterq-dump quit with error code 3

### End Config ###

# open the project specific dataframe, and unroll lists.
def loadProjectData(project_acc,metaDirec = "MetaData" ) : 
    dfpath = metaDirec + "/" + project_acc +"_MetaData.tsv"

    df = pd.read_csv(dfpath, sep= "\t" )
    #### There's probably a better way to do this #####
    # unwrap the Runs, species, taxon 
    dfout = pd.DataFrame()
    for  i in range(df.shape[0]) :
        runs = df.Runs.loc[i].replace("['","").replace("']", "")
        runs = runs.split("', '")

        organism = df.Organism.loc[i].replace("['","").replace("']", "")
        organism = organism.split("', '")
        
        taxid = df.Taxon_ID.loc[i].replace("['","").replace("']", "")
        taxid = taxid.split("', '")
        
        pdat = df.Date.loc[i].replace("['","").replace("']", "")
        pdat = pdat.split("', '")
        
        tmpdf = pd.DataFrame( )
        tmpdf['Run'] = runs
        tmpdf["Submission"] = df.Submission.loc[i]
        tmpdf["Project"] = df.Project.loc[i]
        tmpdf['Experiment'] = df.Experiment.loc[i]
        tmpdf["Organism"] = organism
        tmpdf["Taxon_ID"] = taxid
        tmpdf["Date"] = pdat
        tmpdf["Method"] = df.Method.loc[i]
        
        dfout = pd.concat([dfout , tmpdf],ignore_index=True)

    return (dfout)

# given 10x runs, figure out what version of 10x was used by length of spots
def getReadFilesIn_10x(run , project_acc, fastqstoredir) :


    # look at the first spot - 
    # can edit to avoid writing to disk.
    tmpFile = fastqstoredir+"/"+run+"__temp__.tsv"
    os.system("fastq-dump -X 1 -Z --split-spot "+ run +" > " + tmpFile) 
    with( open(tmpFile ,'r') ) as spot : 
        dat = spot.readlines()
        N = len(dat) //4    # number of fastq files that will be pulled.

    os.remove(tmpFile)

    lengths = {}
    it=1
    for line in dat[0::4] : # look at every 4 lines for the length of the read
        lengths[run +"_"+str(it) +".fastq"] = line.split("length=")[-1].replace("\n","")
        it += 1


    # which read is the longest? 
    l =  [  int(l) for l in lengths.values()  ] 
    ind = l.index(max(l))
    read_bio = list(lengths.keys())[ind]

    # which read contains the barcode + UMI?
    # l =  list(lengths.values() )
    tech = "unknown"
    for i in range(len(l)) :
        n = int(l[i]) 
                
        if n == 24 : 
            tech = "10xv1"
            ind2=i
        elif n == 26 :
            ind2=i
            tech = "10xv2"
        elif n == 28 :
            ind2=i
            tech = "10xv3"
            
    # is tech still unknown? Skip it. 
    if tech =="unknown" :  # doesn't match.. guess the tech by length of cb+umi
        print("... Length of UMI+CB does not match the standard for 10x v1, v2, or v3 for run: "+run+".")
        return (None, None)

        # print("... Run: "  + run +" is not 10x - unable to get technology.")
    else : # tech was found to be a specific 10x version
        read_tech = list(lengths.keys())[ind2]
        readsIn = fastqstoredir+"/"+project_acc+"/"+read_bio + " " +fastqstoredir+"/"+project_acc+"/"+ read_tech

    if len(set(lengths.values())) < 2 : # lengths of reads are identical, unable to distinguish cDNA / CB+UMI
        print("... Lengths of all reads are identical. Unable to distinguish cDNA from CB+UMI for run: "+run+".")
        return (None, None)


    return(readsIn , tech)

# given the 10x df, download the fastq files. Create new file to indicate fastqs are all downloaded.
def getFASTQfiles_10x( df ,fastqstoredir="FASTQ" ,nCore = 30,prefetch =True) :

    # 3210.6243612766266 False - but didn't download any fastq...
    # 5150.785088777542 True - downloaded properly
    # 8779.405247211456 False 
    
    # ti = time.time()
    try :   # make the directory if not already there.
        os.makedirs(fastqstoredir)
    except OSError as err: 
        pass


    # parallelize this. partition nCores / df.shape[0] cores for each run
    # nCoreDump =nCore // df.shape[0]
    nCoreDump = nCore    # temp - pre-parallelized
    AllReadsIn = []
    allTechs= []

    for i in range(df.shape[0]) :
        SRR = df.Run.loc[i]
        SRP = df.Project.loc[i]
        fastqdirec = fastqstoredir+"/"+SRP

        readsIn, tech  = getReadFilesIn_10x(run=SRR,project_acc = SRP,fastqstoredir=fastqstoredir)
        AllReadsIn.append(readsIn)
        allTechs.append(tech)

        isThere = [ os.path.exists(f) for f in readsIn.split(" ")]



        if readsIn != None  and not any (isThere) :
            if prefetch :
                fastqprefix = fastqdirec +"/"+SRR
                print( "... Prefetch SRA file for ",SRR)
                os.system("prefetch --max-size 30000000 -O "+fastqdirec+ " "+ SRR  ) 


                print("... Dumping FASTQ files for "+ SRR)
                os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCoreDump) + " -O "+fastqdirec+ " "+ fastqprefix +".sra" )
            else : 
                print("... Dumping FASTQ files for "+ SRR)
                os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCoreDump) + " -O "+fastqdirec+ " "+ SRR )

        # clean up sra files.
        files = glob.glob(fastqdirec +'/*.sra')

        for f in files:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))

    # tf = time.time()
    # print(tf - ti , prefetch)
    return (AllReadsIn, allTechs)



def buildSmartSeqManifest(project_acc="SRP288492" ,fastqdir="FASTQ" ,metaDirec ="MetaData") : 

    fastqdir = fastqdir + "/" + project_acc
    # list of all fastq files found in the directory
    fastq1files = glob.glob(fastqdir + "/*_1.fastq" )
    fastq2files = glob.glob(fastqdir + "/*_2.fastq" )
    fastqfiles = glob.glob(fastqdir + "/*.fastq" )


    # get the files that dont have a _1 or _2
    fastqfiles = np.setdiff1d(np.setdiff1d(fastqfiles,fastq2files),fastq1files)
    runs = [file.split("_")[0].replace(".fastq","").replace(fastqdir+"/","") for file in fastqfiles ]

    df = pd.DataFrame({"read1":fastqfiles , "read2" :"-", "cellID" : runs} )

    if fastq1files != None  and  fastq2files != None  :
        # contains paired reads. 
        pairedruns = [file.split("_")[0].replace(".fastq","").replace(fastqdir+"/","") for file in fastq1files ]
        df2 = pd.DataFrame({"read1":fastq1files , "read2" :fastq2files, "cellID" : pairedruns} )
        df = pd.concat([df,df2])

        fname = metaDirec + "/"+project_acc+"_SS_manifest.tsv"
        df.to_csv( fname, sep="\t", index=False, header=False) 
    return (df)

# given the SS df, download the fastq files. Create new file to indicate fastqs are all downloaded.
def getReadFilesIn_SS(df, fastqstoredir="FASTQ",metaDirec = "MetaData",nCore = 30) :
    # download the first spot of all reads. 
        
    try :   # make the directory if not already there.
        os.makedirs(fastqstoredir)
    except OSError as err:
        pass

    runList = df.Run.values
    
    # this will be unique 
    project_acc = df.Project.values[0]

    fastqprefix = "/".join([fastqstoredir , project_acc] )

    # batch this 
    start = 0
    batchSize = 100
    n = len(runList)
    while start < n :
        print(start)
        endpoint = min(start+batchSize, n)
        rl = " ".join(runList[start:endpoint] )
        rl2 = [fastqstoredir+"/"+ s+".sra" for s in runList[start:endpoint] ]

        # # this runs sequentially. Is this necessary? Significantly faster with this
        os.system("prefetch --max-size 30000000 -O "+fastqprefix + " "+ rl  ) 
            
        # only get the SRA files that downloaded properly
        sraFiles = glob.glob(fastqprefix+"/*.sra" )
        sraFiles =  list(set(sraFiles)& set(rl2) )
        sraFiles =" ".join(sraFiles)

        # which runs don't have a fastq?
        # fastqFiles =  glob.glob(fastqstoredir+"/*fastq" )
        os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCore) + " -O "+fastqprefix+  " "+ rl )

        start = min(start + batchSize , n)
        # clean up the sra files. 
        # os.remove(fastqstoredir + "*.sra")

        # clean up sra files.
        files = glob.glob(fastqprefix +'/*.sra')

        for f in files:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))

    # project is done - build the manifest
    df = buildSmartSeqManifest(project_acc ,fastqstoredir ,metaDirec)
    return(df )

#####################################


# this just downloads the FASTQ files 
def main( project_acc ="SRP288492", species ="mouse",metaDirec ="MetaData",fastqstoredir ="FASTQ", nCore=30) :
    
    df = loadProjectData(project_acc, metaDirec )

    # split the data into 10x portion and SS portion. 
    dfByMethod ={ x:y for x,y in  df.groupby("Method") }
    
    key = list(dfByMethod.keys())
    if "isSome10x" in key :
        df10x = dfByMethod["isSome10x"]    
        AllReadsIn, allTechs = getFASTQfiles_10x(df10x, fastqstoredir,nCore ,True)
        # AllReadsIn2, allTechs2 = getFASTQfiles_10x(df10x, fastqstoredir ,nCore ,True)
        # Store these to a file indicating the fastq files for the project are all downloaded. 

        #which runs failed? 
        # df10x.loc[ np.where(allTechs == None) ,:]
        filename = metaDirec + "/" + project_acc + "_10x_MetaData.tsv" 
        df10x["readsIn"]= AllReadsIn
        df10x["Method"] =  allTechs
        df10x.to_csv(filename , sep ="\t" ,index=False)

        # Separately run star for all of these readsIn


    if "isSS" in key :
        dfSS = dfByMethod["isSS"]    
        manifestPath =getReadFilesIn_SS(dfSS, fastqstoredir ,nCore)
    
        filename = metaDirec + "/Projects/" + species + "/" + project_acc + "_SS_MetaData.tsv" 
        # df10x["readsIn"]= manifestPath
        dfSS["Method"] =  "SMART-Seq"

        dfSS.to_csv(filename , sep ="\t" ,index=False)
        maniname = metaDirec + "/Projects/" + species + "/" + project_acc + "_SS_STAR_manifest.tsv" 
        manifestPath.to_csv(maniname,  sep ="\t" ,index=False)



main( project_acc , species ,metaDirec ,fastqStore , nCore)




