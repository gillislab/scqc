# get FASTQ files
# import pickle
import contextlib
import os 
import subprocess
import sys
import h5py
import pandas as pd
import glob
from contextlib import contextmanager
# part 3: Get the fastq files for the specified project ID, verify technology and align with STAR

# use:
# python getFastq2STAR.py "SRX7719269" "mouse" "SRA_mouse_scRNA.pkl" "SupplementData" "FASTQ" "STARout" "30"

project_acc = sys.argv[1]
# species = sys.argv[2]
filepath = sys.argv[2]
suppDirec = sys.argv[3]
fastqStore = sys.argv[4]
STARout = sys.argv[5]
nCore=sys.argv[6]

project_acc = 'SRP288492'


fastqStore="FASTQ"
nCore = 30      # max cores allows
filepath="MetaData/SRA_MetaData_mouse_all_Current.hdf5"
suppDirec = "SupplementData"
# species = "mouse"
STARout = "STARout"



tax2spec = {
    "10090":"mouse",
    "9606":"human"
}

@contextmanager
def suppress_stderr():
    with open(os.devnull, "w") as devnull:
        old_stderr = sys.stderr
        sys.stderr = devnull
        try:  
            yield
        finally:
            sys.stderr = old_stderr

def loadProjectData(project_acc,filepath ) : 

    allRows =[]
    with h5py.File(filepath ,'r') as f : 
        for srx in f[project_acc].keys() : 
            ky = "/".join([project_acc,srx])
            runs = f[ky]['Runs'][()] 
            method = f[ky]['Method'][()] 
            taxid = f[ky]['Taxon_ID'][()] 

            row = [project_acc, srx, runs , method, taxid ]
            allRows.append(row)

    df = pd.DataFrame(allRows,columns = ["Project","Experiment","Runs","Method","Taxon_ID" ] )


    return (df)

def getReadFilesIn_10x(run,  exp_dic, fastqstoredir ) :

    # look at the first spot - 
    os.system("fastq-dump -X 1 -Z --split-spot "+ run +" > __temp__.tsv") 
    with( open("__temp__.tsv" ,'r') ) as spot : 
        dat = spot.readlines()
        N = len(dat)

    os.remove("__temp__.tsv")

    lengths = {}
    it=1
    for line in dat[0::4] : # look at every other line for the length of the read
        lengths[run +"_"+str(it) +".fastq"] = line.split("length=")[-1].replace("\n","")
        it += 1


    # which read is the longest? 
    l =  list(lengths.values() )
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

        

    # does the length of the read match what we expect?
    if exp_dic["Method"] ==tech : # length of cb + umi matches the identified technology.
        read_tech = list(lengths.keys())[ind2]
        readsIn = fastqstoredir+"/"+read_bio + " " +fastqstoredir+"/"+ read_tech

    elif tech =="unknown" :  # doesn't match.. guess the tech by length of cb+umi
        readsIn = " "
        print("... Run: "  + run +" is not 10x - unable to get technology.")
    else :
        print('... Length of UMI+CB does not match the default for '+exp_dic["Method"]+ " using "+ tech+ " instead")
        read_tech = list(lengths.keys())[ind2]
        readsIn = fastqstoredir+"/"+read_bio + " " +fastqstoredir+"/"+ read_tech

    if read_tech == read_bio :
        readsIn = " "
        print("... No biological reads found! ")

    return(readsIn , tech)


df = loadProjectData(project_acc,filepath )
# split the data into 10x portion and SS portion. 
ssDat = df.loc[df.Method == "isSS",:]
tenxDat =df.loc[df.Method.str.contains("10x"),:]

runList = ssDat.Runs.values
fastqstoredir = "/".join([fastqStore, project_acc])

def getReadFilesIn_SS(runList, fastqstoredir,nCore = 5) :
    # download the first spot of all reads. 
        
    try :   # make the directory if not already there.
        os.makedirs(fastqstoredir)
    except OSError as err: 
        pass


    # parallelize this - very slow as is. ~ 2 sec / run
    # do we need to do this? Why not just download the fastq files?
    # maniList = []
    # for run in runList :
    #     # tmp = os.popen("fastq-dump -X 1 -Z --split-spot "+ run).read().split("\n")
    #     # print(run)
    #     tmp = subprocess.Popen("fastq-dump -X 1 -Z --split-spot "+ run, shell=True ,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0].decode().split("\n")
    #     tmp = tmp[0::4]
    #     for i in range(len(tmp)) :
    #         a= tmp[0].split("length=")[-1].isdigit()           
    #         b= tmp[1].split("length=")[-1].isdigit()           

    #     if a & b :  # both reads file found      
    #         row = [fastqstoredir+"/"+run+"_1.fastq" , fastqstoredir+"/"+run+"_2.fastq", run ] 
    #         maniList.append(row)
    #     elif a and not b:    # only one file found
    #         row = [fastqstoredir+"/"+run+".fastq" , "-", run ] 
    #         maniList.append(row)
    #     elif not a and not b : # no spot found
    #         print("... ... No spots found for " + run)

    # manifest = pd.DataFrame(maniList, columns= ["read_1" , "read_2","cell_id" ]  )
    
    
    
    # batch this 
    start = 0
    batchSize = 100
    n = len(runList)
    while start < n :
        print(start)
        endpoint = min(start+batchSize, n)
        rl = " ".join(runList[start:endpoint] )
        rl2 = [fastqstoredir+"/"+ s+".sra" for s in runList[start:endpoint] ]

        # this runs sequentially.
        os.system("prefetch --max-size 30000000 -O "+fastqstoredir + " "+ rl  ) 
            
        # only get the 
        sraFiles = glob.glob(fastqstoredir+"/*.sra" )
        sraFiles =  list(set(sraFiles)& set(rl2) )
        sraFiles=" ".join(sraFiles)

        # which runs don't have a fastq?
        # fastqFiles =  glob.glob(fastqstoredir+"/*fastq" )
        os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCore) + " -O "+fastqstoredir+  " "+ sraFiles )
        start = min(start + batchSize , n)
        # clean up the sra files. 
        # os.remove(fastqstoredir + "*.sra")

#####################################


    


    

def getFASTQfile( runlist ,outdir) :

    pass




try :
    fastqstoredir=fastqStore+"/"+project_acc
    os.makedirs(fastqstoredir)
except OSError as error: 
    pass

try :
    staroutdir=STARout+"/"+project_acc
    os.makedirs(staroutdir)
except OSError as error: 
    pass



if ssDat.shape[0]  > 0 : # at least one experiment found containing smartseq data
    # download all the files for all runs 

    pass

if tenxDat.shape[0]  > 0  : # at least one experiment found containing 10x data.
    pass






genomeDir = suppDirec + "/genomes/" + species + "/STARindices"

if  not os.path.isdir(genomeDir)  :
    raise NameError("Genome Directory does not exist: "+ genomeDir) 

















SRXdict = loadSRA_MetaData(filepath)
exps = SRXdict.keys()

sra = SRXdict[exp]["Submission"]
runs = list(SRXdict[exp]["Runs"])

nCoreDump = int(nCore / len(runs))



params = {}
params["10xv1"] = {        
    "solo_type":"CB_UMI_Simple",
    "whiteListPath" : suppDirec+"/whiteLists/737K-april-2014_rc.txt",
    "CB_length":"14",
    "UMI_start":"15",
    "UMI_length":"10" 
}

params["10xv2"] = {        
    "solo_type":"CB_UMI_Simple",
    "whiteListPath" : suppDirec+"/whiteLists/737K-august-2016.txt",
    "CB_length":"16",
    "UMI_start":"17",
    "UMI_length":"10" 
}

params["10xv3"] = {        
    "solo_type":"CB_UMI_Simple",
    "whiteListPath" : suppDirec+"/whiteLists/3M-february-2018.txt",
    "CB_length":"16",
    "UMI_start":"17",
    "UMI_length":"12" 
}

params["SS"] ={}

    # if ( tech == "SmartSeq" ){
    #     solo_type="SmartSeq"
    #     path2manifest = paste0(manifests)
    # } 


# look for  10x runs - do these first.
SSruns = [] # if we find smart seq runs, append them to a list and build the manifest after.
for run in runs :   # parallelize this. 
    # make sure the Method field is there and is allowed
    flag=0 # unsure if we should run

    try :  
        SRXdict[exp]["Method"]
        if SRXdict[exp]["Method"]  in ['10xv1','10xv2','10xv3',"Some10x"] : 
            flag = 1  
        elif SRXdict[exp]["Method"]  == "SS" : 
            flag = 2
               
    except KeyError :
        flag = 0


    if flag == 1 :  # 10x chromium
        
        readsIn, tech  = getReadFilesIn_10x(run,  SRXdict[exp] , fastqstoredir)

        print( "... Prefetch SRA file for ",run)
        os.system("prefetch --max-size 30000000 -O "+fastqStore+"/"+sra+"/"+exp + " "+ run  ) 


        print("... Dumping FASTQ files for "+ run)
        os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCoreDump) + " -O "+fastqStore+ "/"+sra+"/"+exp+ " "+ fastqStore+ "/"+sra+"/"+exp+"/"+run+".sra" )

        if tech != "unknown" :
            param = params[tech ]

            # RUN STAR USING THE FASTQ FILES
            os.system( (
                "STAR" + 
                " --runMode alignReads" +
                " --runThreadN "+ str(nCore)+
                " --genomeDir "+ genomeDir+
                " --outFileNamePrefix "+ staroutdir+ "/" + run  +
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
        else : 
            print("Technology is unknown - skipping "+ run)



        # else :          # unsure - get the BAM files
            # os.system("sam-dump --unaligned --output-file "+fastqStore+"/"+sra+"/"+run+".sam " + fastqStore+ "/"+sra+"/"+run+".sra"  )

    
    elif flag == 2 : # smart seq 
        
        print( "... Prefetch SRA file for ",run)
        os.system("prefetch --max-size 30000000 -O "+fastqStore+"/"+sra+"/"+exp + " "+ run  ) 


        print("... Dumping FASTQ files for "+ run)
        os.system("fasterq-dump --split-files --include-technical --threads "+ str(nCoreDump) + " -O "+fastqStore+ "/"+sra+"/"+exp+ " "+ fastqStore+ "/"+sra+"/"+exp+"/"+run+".sra" )

        pass


