# get all UIDs for mouse RNA seq experiments
import datetime as dt
import json
import os
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

getAll = False          # Boolean 
getRecent = True        # Boolean - ignored if getAll == True
lastNdays = 5           # only used if getRecent == True
metaDir = "MetaData"    # 
### End Config ###

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

    # can we get all of the UIDs found in the search? if not, loop through the rest
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

def fetchData(UIDList,start = 0, batchSize= 200, outpath="MetaData") :
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

        endpoint = min(n, start+batchSize)
        
        try : 
            uid_str= ",".join([str(UID) for UID in UIDList[start:endpoint ] ])
            http = urllib3.PoolManager()
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="+uid_str
            r = http.request("POST", url)
            r = r.data.decode()
            root = ET.fromstring(r)
            # write to csv
            pd.DataFrame( UIDList[start:endpoint ]).to_csv(outpath+"/UIDs_fetched.csv",mode="a", header=False, index=False)
            flag = 0
        except:
            print('... ... ...Broken pipe at ' + str(start) + " : ", str(endpoint))
            # brokenPipes[start] =  endpoint

            # write to csv
            pd.DataFrame( UIDList[start:endpoint ]).to_csv(outpath+"/UIDs_Not_fetched.csv",mode="w", header=False, index=False)

            flag = -1 # fail flag - don't continue

        
        if flag == 0 :
            for exp in root.iter("EXPERIMENT_PACKAGE"):

                # lcps = []
                for lcp in exp[0].iter("LIBRARY_CONSTRUCTION_PROTOCOL") :
                    lcp = lcp.text

                
                SRXs = exp[0].get('accession')
                SRAs = exp[1].get('accession')
                SRPs = exp[3].get('accession')
                # title = exp[3][1][0].text
                abstract =exp[3][1][2].text
                SRRs = []
                date=[]
                taxon=[]
                orgsm =[]
                for study in exp[3].iter("STUDY_TITLE") :
                    title = study.text 
                for study in exp[3].iter("STUDY_ABSTRACT") :
                    abstract = study.text 

                for run in exp.iter("RUN") : 
                    SRRs.append(run.attrib['accession'])
                    date.append(run.attrib['published'])
                    for mem in run.iter("Member") : 
                        taxon.append(mem.attrib['tax_id'])
                        orgsm.append(mem.attrib['organism'])
                

                    
                row = [SRPs, SRXs, SRAs, SRRs, date, taxon,orgsm, lcp, title,abstract]
                allRows.append(row)
                # {0:"Project",1: "Experiment",2: "Submission", 3:"Runs",4: "Taxon_ID", 5:"LCP" }

                # experiments[SRXs] = {"LibConsProt": lib_cons_prots , "Runs":SRRs ,"Submission":SRAs,"Project":SRPs }


        start = start + batchSize
        

        if  (it % 100  == 0)  :
            # try : 
            df = pd.DataFrame(allRows , columns = ["Project","Experiment","Submission", "Runs","Date","Taxon_ID", "Organism","LCP","Title","Abstract"  ])
            df = df.fillna(value = "")  # fill None with empty strings. 


            if os.path.exists(outpath+"/allMetaData.tsv") : 
                #path exists - remove header then append
                keepHeader = False
            else :  # path doesn't exist, keep header
                keepHeader=True
            
            df.to_csv(outpath+"/allMetaData.tsv" ,sep="\t", mode = 'a', index=False,header=keepHeader)

            print( "... ..."+ str(round(it/Nbatch*100,2))  +"% done" )
            saveAsFiles(df,outpath)
            # saveAsHDF5(df ,outpath)
            allRows =[]



        it +=1


        


    df = pd.DataFrame(allRows , columns = ["Project","Experiment","Submission", "Runs","Date","Taxon_ID", "Organism","LCP","Title","Abstract"  ])
    df = df.fillna(value = "")  # fill None with empty strings. 

    if os.path.exists(outpath+"/allMetaData.tsv") : 
        #path exists - remove header then append
        keepHeader = False
    else :  # path doesn't exist, keep header
        keepHeader=True
    
    df.to_csv(outpath+"/allMetaData.tsv" ,sep="\t", mode = 'a', index=False,header=keepHeader)
    saveAsFiles(df, outpath)

    return( )

# not used
def saveAsHDF5(df,outpath = "MetaData",add2existing=True  ):
    # if h5File == None : 
    #     h5File = "SRA_data.hdf5"
    if add2existing : 
        file = "current"
    else : 
        file = str(dt.date.today())


    h5File = outpath + "/SRA_data_"+file+".hdf5"
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

def saveAsFiles(df, outpath = "MetaData") : 
    df['Method'] = "Unknown"
    df["Status"] = "UIDfecthed"
    # split the dataframe by project
    by_usrp = [y for x, y in df.groupby('Project', as_index=False)]
    for i in range(len(by_usrp)) : 
        # does the file exist? 
        srp = by_usrp[i].Project.unique()[0]
        if os.path.exists(outpath+"/"+srp+"/MetaData.tsv") : 
            # create directory
            #path exists - remove header then append
            keepHeader = False
        else :  # path doesn't exist, keep header and make the directory
            keepHeader=True
            os.makedirs(outpath+"/"+srp)

        by_usrp[i].to_csv(outpath+"/"+srp+"/MetaData.tsv" , sep="\t" , mode= "a" ,index=False, header=keepHeader)
    return ()    

def main(species = "mouse",getAll=False, getRecent = True,lastNdays=5, metaDir = "MetaData" ):

    assert(getAll | getRecent )
    assert(species in ["mouse","human"] )
    assert(type(lastNdays) == int )
    assert(type(metaDir) == str)

    try :
        os.makedirs(metaDir)
    except OSError : 
        pass

    ti=time.time()

    if getAll:
        s = "all"
        statement='\n Getting list of '+s+" " +species+' single cell RNA-seq UIDs...'

    elif getRecent : 
        s = "recent"
        s2 = str(lastNdays)
        statement='\n Getting list of '+s+" " +species+' single cell RNA-seq UIDs in the last '+s2+" days..."

    print(statement)

    batchSize = 200     # max 200 
    retmax = 100000     # max 100000 - how many should we gather at a time?


    UIDList = getUIDlist(species , getAll , getRecent , lastNdays, retmax )

    # which UIDs failed previously? Let's add those in
    if os.path.exists(metaDir +"/UIDs_Not_fetched.csv"):
        UIDsnotFetched = pd.read_csv(metaDir +"/UIDs_Not_fetched.csv" ,header=None)[0].values
        UIDList =  np.append(UIDsnotFetched ,UIDList)

    # which UIDs did we already fetch? - lets not do those again    
    if os.path.exists(metaDir +"/UIDs_fetched.csv" ) :
        # Get the UIDs in UIDList but not in the fetched file
        UIDsDone = pd.read_csv(metaDir +"/UIDs_fetched.csv" ,header=None)[0].values
        # which UIDList are not in UIDsDone?
        UIDList = np.setdiff1d(UIDList,UIDsDone)

    

    if len(UIDList) <= 0 :
        return
    
    
    print('... Converting these to SRA accession IDs and getting MetaData')
    fetchData(UIDList,start = 0, batchSize=batchSize, outpath=metaDir )


    # trying again with the ones that we missed, decreasing the batch size with each iteration
    print('... Going back to the ones that we missed.')

    if os.path.exists(metaDir +"/UIDs_Not_fetched.csv"):
        UIDs = pd.read_csv(metaDir +"/UIDs_Not_fetched.csv" ,header=None)[0].values
        lold = len(UIDs)
        while lold > 0 : 
            if batchSize <= 1 :
                break
            batchSize = round(batchSize / 2)
            fetchData(UIDs,start = 0, batchSize=batchSize, outpath=metaDir )
            UIDs = pd.read_csv(metaDir +"/UIDs_Not_fetched.csv" ,header=None)[0].values
            lnew = len(UIDs)

            if lnew == lold :  # All UIDs are the same as before.. unable to get more
                break
            
    tf = time.time()

    print(tf-ti)






main(species, getAll, getRecent, lastNdays, metaDir)
