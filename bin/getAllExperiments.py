# get all UIDs for mouse RNA seq experiments
import datetime as dt
import json
import os
import pickle
import time
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd
import urllib3

# part 1 - get the SRA meta data. Either for all datasets or for last few days.

### Config ####
# getUID options
species = "mouse"   # mouse or human for now           
batchSize = 200     # max 200 
retmax = 100000     # max 100000 - how many should we gather at a time?
getAll = True       # Boolean 
getRecent = True    # Boolean - ignored if getAll == True
lastNdays = 5       # only used if getRecent == True
metaDir = "MetaData"
saveIt=True

### End Config ###

try :
    os.makedirs(metaDir)
except OSError as error: 
    pass


if getAll :
    st = "all"
elif  getRecent :
    st = "Last_"+str(lastNdays)+"days"

today = str(dt.date.today())

SRAmd_file = metaDir+"/SRA_MetaData_"+species+"_"+st +"_"+today+".hdf5"
noLCP = metaDir+"/SRA_MetaData_"+species+"_"+st +"_"+today+"_BrokenPipes.pkl"
currentFile =metaDir+"/SRA_MetaData_"+species+"_"+st +"_current.hdf5"



spec2sci = {
    "mouse":"mus+musculus",
    "human":"homo+sapiens"
}

spec2taxon = {
    "mouse":"10090",
    "human":"9606"
}


def getUIDlist ( species = "mouse", getAll = True, getRecent = False, lastNdays =5, retmax =100000  ):
    # defaults to getting RNA-Seq experiments with key word "single cell"
    retstart=0
    spec2sci = {
        "mouse":"mus+musculus",
        "human":"homo+sapiens"
    }

    http = urllib3.PoolManager()

    if  getAll :
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%22rna+seq%22[Strategy]+%22"+spec2sci[species]+ "%22[Organism]+%22single+cell%22[Text Word]&retstart="+str(retstart)+"&retmax="+str(retmax)+"&retmode=json" 
    elif getRecent :
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%22rna+seq%22[Strategy]+%22"+spec2sci[species]+"%22[Organism]+%22single+cell%22[Text Word]&retstart="+str(retstart)+"&retmax="+str(retmax)+"&datetype=pdat&reldate="+str(lastNdays)+"&retmode=json"
    
    r = http.request("GET", url)
    esearch_res = json.loads(r.data.decode('utf-8'))
    nUIDs   = int(esearch_res['esearchresult']['count'] )
    UIDs    = esearch_res['esearchresult']['idlist']

    # get we get all of the UIDs found in the search? if not, loop through the rest
    it=1

    while nUIDs >= (len(UIDs)+3) :
        nLast = len(UIDs) 
        retstart_new = retmax*it + 1
        if  getAll :
            url_new = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%22rna+seq%22[Strategy]+%22"+spec2sci[species]+ "%22[Organism]+%22single+cell%22[Text Word]&retstart="+str(retstart_new)+"&retmax="+str(retmax)+"&retmode=json" 
        elif getRecent :
            url_new = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%22rna+seq%22[Strategy]+%22"+spec2sci[species]+"%22[Organism]+%22single+cell%22[Text Word]&retstart="+str(retstart_new)+"&retmax="+str(retmax)+"&datetype=pdat&reldate="+str(lastNdays)+"&retmode=json"


        rtmp = http.request("GET", url_new)
        esearch_res_tmp = json.loads(rtmp.data.decode('utf-8'))
        # nUIDs_tmp   = int(esearch_res['esearchresult']['count'] )
        UIDs_tmp    = esearch_res_tmp['esearchresult']['idlist']

        UIDs = UIDs + UIDs_tmp
        
        it += 1

        if nLast == len(UIDs) :
            break
        
    print("... We found " +str(len(np.unique(np.asarray(UIDs))))+ " UIDs." )

    return ( np.unique(np.asarray(UIDs)) )

def fetchData(UIDList,start = 0, batchSize= 200, saveIt=True, outfile="MetaData/SRA_MetaData.hdf5",failFile = "MetaData/brokenPipes.pkl") :
    # outfile is ignored if saveIt is False
    # batch it
    n = len(UIDList) 
    # batchSize = 200
    # start = 0
    Nbatch = np.ceil(n/batchSize)
    brokenPipes = {}
    # experiments = {}
    it = 1
    allRows =[]

    while  start <= n :

        # print(start)
        # print(it)
        endpoint = min(n, start+batchSize)
        
        try : 
            uid_str= ",".join([str(UID) for UID in UIDList[start:endpoint ] ])
            http = urllib3.PoolManager()
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="+uid_str
            r = http.request("POST", url)
            r = r.data.decode()
            root = ET.fromstring(r)
            flag = 0
        except:
            print('... ... ...Broken pipe at ' + str(start) + " : ", str(endpoint))
            brokenPipes[start] =  endpoint
            flag = -1 # fail flag - don't continue

        
        if flag == 0 :
            for exp in root.iter("EXPERIMENT_PACKAGE"):

                # lcps = []
                for lcp in exp[0].iter("LIBRARY_CONSTRUCTION_PROTOCOL") :
                    lcp = lcp.text

                
                SRXs = exp[0].get('accession')
                SRAs = exp[1].get('accession')
                SRPs = exp[3].get('accession')

                SRRs = []
                date=[]
                taxon=[]
                orgsm =[]
                for run in exp.iter("RUN") : 
                    SRRs.append(run.attrib['accession'])
                    date.append(run.attrib['published'])
                    for mem in run.iter("Member") : 
                        taxon.append(mem.attrib['tax_id'])
                        orgsm.append(mem.attrib['organism'])
                    
                

                    
                row = [SRPs, SRXs, SRAs, SRRs, date, taxon,orgsm, lcp ]
                allRows.append(row)
                # {0:"Project",1: "Experiment",2: "Submission", 3:"Runs",4: "Taxon_ID", 5:"LCP" }

                # experiments[SRXs] = {"LibConsProt": lib_cons_prots , "Runs":SRRs ,"Submission":SRAs,"Project":SRPs }


        start = start + batchSize
        

        if  (it % 100  == 0)  & saveIt:
            # try : 
            df = pd.DataFrame(allRows , columns = ["Project","Experiment","Submission", "Runs","Date","Taxon_ID", "Organism","LCP"  ])
            df = df.fillna(value = "")  # fill None with empty strings. 

            print( "... ..."+ str(round(it/Nbatch*100,2))  +"% done" )
            
            saveAsHDF5(df ,outfile)
            allRows =[]

            remainingUIDs = UIDList[endpoint:-1]
            UIDsDone = UIDList[:endpoint]
            pd.DataFrame(remainingUIDs ).to_csv("remainingUIDs.csv",header=False,index=False)
            pd.DataFrame(UIDsDone ).to_csv("UIDsDone.csv",header=False,index=False)


            # f = open(outfile ,"wb")
            # pickle.dump(df,f)
            # f.close()

            f2 = open(failFile,"wb")
            pickle.dump(brokenPipes,f2)
            f2.close()
            


        it +=1


        
    if saveIt :  

        df = pd.DataFrame(allRows , columns = ["Project","Experiment","Submission", "Runs","Date","Taxon_ID", "Organism","LCP"  ])
        df = df.fillna(value = "")  # fill None with empty strings. 

        saveAsHDF5(df ,outfile)

        print("Done. Saving to ",outfile)
        # f = open(outfile ,"wb")
        # pickle.dump(df,f)
        # f.close()

        remainingUIDs = UIDList[endpoint:-1]
        UIDsDone = UIDList[:endpoint]
        pd.DataFrame(remainingUIDs ).to_csv("remainingUIDs.csv",header=False,index=False)
        pd.DataFrame(UIDsDone ).to_csv("UIDsDone.csv",header=False,index=False)


        f2 = open(failFile,"wb")
        pickle.dump(brokenPipes,f2)
        f2.close()
        
        
    else :
        print("Done. Not saving to a file.")


    
    return( df , brokenPipes)

def saveAsHDF5(df, h5File = None  ):
    if h5File == None : 
        h5File = "SRA_data.hdf5"


    projs = np.unique(df.Project.values)

    with h5py.File(h5File,"a")  as f : 

        for proj in projs :    # project IDs
            if proj in f.keys() :
                srp = f[proj]
            else :
                srp = f.create_group(proj)

            # for each project, look at the unique experiments.
            dftmp = df.loc[df.Project == proj,: ]


            for exp in np.unique(dftmp.Experiment.values):
                tmp =dftmp.loc[dftmp.Experiment == exp ,:]
                if exp in srp.keys() :  
                    srxDat  = srp[exp]
                else :
                    srxDat = srp.create_group(exp)

                try : 
                    srxDat.create_dataset("Submission", data = tmp.Submission.values[0] )
                    srxDat.create_dataset("LCP" , data = tmp.LCP.values[0] )
                    srxDat.create_dataset("Runs" , data = ",".join(tmp.Runs.values[0]) )
                    srxDat.create_dataset("Date" , data = ",".join(tmp.Date.values[0]) )
                    srxDat.create_dataset("Taxon_ID" , data = ",".join(tmp.Taxon_ID.values[0]) )
                    srxDat.create_dataset("Organism" , data = ",".join(tmp.Organism.values[0]) )
                except :
                    print("skipping " +proj + "/"+exp  +'/. already exists')



                # # tmp.Runs.values
                # for  i in range(len(tmp.Runs.values[0])) : 
                #     # if tmp.Runs.values[0][i] in 
                #     srrDat= srxDat.create_group( tmp.Runs.values[0][i] )
                #     srrDat.create_dataset("Date" , data = tmp.Date.values[0][i])
                #     srrDat.create_dataset("Taxon_ID" , data = tmp.Taxon_ID.values[0][i])

            
            # to access data. 
            # f['SRP285322/SRX9187466/SRR12708542/Date/'][()]
            # f['SRP285322/SRX9187466/LCP'][()]


    pass


def loadPickle(filepath): 
    
    with (open(filepath, "rb")) as openfile:
        while True:
            try:
                a=pickle.load(openfile)
                return (a)

            except EOFError:
                pass





ti=time.time()
if getAll:
    s = "all"
    statement='Getting list of '+s+" " +species+' single cell RNA-seq UIDs...'

elif getRecent : 
    s = "recent"
    s2 = str(lastNdays)
    statement='Getting list of '+s+" " +species+' single cell RNA-seq UIDs in the last '+s2+" days..."

print(" ")
print(statement)


UIDList = getUIDlist(species , getAll , getRecent , lastNdays, retmax )
print('... Converting these to SRA accession IDs and getting MetaData')
df , brokenPipes = fetchData(UIDList,start = 0, batchSize=batchSize, saveIt=True,outfile = SRAmd_file,failFile=noLCP)

# start = 5300000
tf=time.time()

print("Finished in " +str((tf -ti) / 60) +" Min")




# continue with the ones that broke. 



# bp = loadSRA_MetaData("MetaData/SRA_MetaData_mouse_all_2021-05-12_BrokenPipes.pkl2")
starts = list(brokenPipes.keys())
ends =  list(brokenPipes.values())


UIDs = []
for i in range(len(starts)) : 
    UIDs.extend(list(UIDList.loc[starts[i]:ends[i],0] ))

while len(UIDs) > 0 : 
    batchSize = round(batchSize / 2)
    df , brokenPipes = fetchData(UIDs,start = 0, batchSize=batchSize, saveIt=True,outfile = SRAmd_file,failFile=noLCP+"2")
    
    starts = list(brokenPipes.keys())
    ends =  list(brokenPipes.values())

    UIDList = UIDs.copy()
    UIDs=[]
    for i in range(len(starts)) : 
        UIDs.extend(UIDList[starts[i]:ends[i]] )

    if len(UIDs) == len(UIDList) :  # All UIDs are the same as before.. unable to get more
        break
    
    if batchSize == 1 :
        break

# Get the corresponding technology iif able
# download the fastq files if paired. 


