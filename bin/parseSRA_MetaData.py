import os

import h5py
import pandas as pd
# from numpy.lib.npyio import save

# import pickle


h5path ="MetaData/SRA_MetaData_mouse_all_Current.hdf5"
saveIt = True

# import sys

# sys.argv

# part 2: predict the technology used from the LCP/ 
# only needs to be done once. 
# this script takes the results from 'getAllExperiments.py' and  builds a dataframe 




def hdf5_to_df(h5path, saveIt = True) : 
    # with h5py.File(h5path,'r') as f : 

    tsvpath= h5path.replace("hdf5","tsv")

    if os.path.isfile(tsvpath) :    # already done - load it
        df = pd.read_csv(tsvpath,sep="\t")
        print("This hdf5 file has already been converted to a dataframe. Loading it instead.")
    else : 
        with  h5py.File(h5path,'r')  as f :
            projs = f.keys()
            n = len(projs)
            it = 1
            allRows = []
            for proj in projs : 
                exps = f[proj].keys()
                # print(proj)
                for exp in exps : 
                    srx = "/".join([proj,exp])

                    fields = f[srx].keys() 

                    row = [proj, exp ]
                    for field in fields  : 
                        row.append(f[srx][field][()])

                    allRows.append(row)

                if it % 100 == 0 : 
                    print("... " + str(round(it/n*100,2)) +"% done." )
                it = it +1 

            df = pd.DataFrame(allRows, columns = ["Project","Experiment" ] + list(fields))

            if saveIt :
                df.to_csv(tsvpath,sep= "\t", index=False)

    return (df)
# conversion is pretty slow. lol            

def predictTechFromLCP(df) :
    # loop through dictionary, search for keywords in 
    keywords = {
        "is10x" :"10x Genomics|chromium|10X protocol|Chrominum|10X 3' gene|10X Single|10x 3'",
        "v3" : "v3 chemistry|v3 reagent|V3 protocol|CG000206|Single Cell 3' Reagent Kit v3|10X V3|1000078",
        "v2" : "v2 chemistry|v2 reagent|V2 protocol|P/N 120230|Single Cell 3' v2|Reagent Kits v2|10X V2",
        "v1" : "Kit v1|PN-120233|10X V1",
        "SS" : "Smart-Seq|SmartSeq|Picelli|SMART Seq",
        "smarter":"SMARTer",
        "dropseq":"Cell 161, 1202-1214|Macosko|dropseq|drop-seq",
        "celseq":"CEL-Seq2|Muraro|Cell Syst 3, 385|Celseq2|Celseq1|Celseq|Cel-seq",
        "sortseq":"Sort-seq|Sortseq|Sort seq",
        "seqwell" : "Seq-Well|seqwell",
        "biorad": "Bio-Rad|ddSeq",
        "indrops":"inDrop|Klein|Zilionis",
        "marsseq2":"MARS-seq|MARSseq|Jaitin et al|jaitin",
        "tang":"Tang",
        # "TruSeq":"TruSeq",
        "SPLiTseq": "SPLiT-seq",
        "Microwellseq":"Microwell-seq"
    }


    uLCP = pd.DataFrame({"LCP" : df.LCP.unique() })
    
    # search for the keywords
    for i in range(len(keywords)) :
        key= list(keywords)[i]
        kw = keywords[key] 
        uLCP[key] = uLCP.LCP.str.lower().str.contains(kw.lower())

    # only keep the ones where a keyword was found
    # ["Experiment","Submission","Runs","Project"]
    # cols = df.columns.values

    tmpDF = uLCP.loc[:,list(keywords.keys())]
    tmpDF = tmpDF.fillna(False)



    unknownTechs = uLCP.loc[tmpDF.sum(axis=1) ==0,"LCP"]
    tmpDF["isSome10x"]  = tmpDF.loc[:,"is10x"]


    tmpDF["isSS"]    = tmpDF.loc[:,"SS"] & ~(tmpDF.loc[:,"is10x"].astype('bool'))
    tmpDF["isMultiple"] = tmpDF.iloc[:,4:-1].sum(axis=1) +tmpDF.iloc[:,0]  > 1 

    tmpDF["is10xv3"] = tmpDF.loc[:,"v3"]
    tmpDF["is10xv2"] = tmpDF.loc[:,"v2"]
    tmpDF["is10xv1"] = tmpDF.loc[:,"v1"]

    tmpDF["Method"] = "Other"

    for tech in ["isSome10x","is10xv1","is10xv2","is10xv3","isSS" ] :
        tmpDF.loc[tmpDF.loc[:,tech],"Method"]  = tech

    
    uLCP['Method'] = tmpDF.Method
    uLCP = uLCP.loc[:,["LCP","Method"]]
    df = df.merge(uLCP , on = "LCP")

    

    return(df, unknownTechs)    



df = hdf5_to_df(h5path,saveIt )
# did we already do this?
try :
    df.Method
    print("Already has a method field")    
except : 
    df , unknownLCP = predictTechFromLCP(df)




if saveIt : 
    # save the df.
    tsvpath= h5path.replace("hdf5","tsv")
    df.to_csv(tsvpath,sep= "\t", index=False) 

    try: 
        unknownLCPFile = tsvpath.replace(".tsv","_unknownTech.tsv")
        unknownLCP.to_csv(unknownLCPFile,sep= "\t", index=False) 
    except: 
        pass

    # save the method in the hdf5 file

    with h5py.File(h5path,"a")  as f : 
        for i in range(df.shape[0]) : 
            ky = "/".join(df.loc[i,["Project","Experiment"] ].values)+"/"
            try: 
                f[ky].create_dataset("Method", data= df.loc[i,"Method"] )
                print('creating ds for '+ ky)
            except OSError as err: 
                # print(str(err) + " '" + ky +"'")
                pass
