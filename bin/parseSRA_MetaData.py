import os

import h5py
import pandas as pd
import glob
# from numpy.lib.npyio import save

# load the full df. 


### config ###
metaDirec = "MetaData"
species = "mouse"
getAll = False
getRecent = True

# part 2: predict the technology used from the LCP/ 
# only needs to be done once. 
# this script takes the results from 'getAllExperiments.py' and  builds a dataframe 




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


    # tmpDF["isMultiple"] = tmpDF.iloc[:,4:-1].sum(axis=1)  > 1 
    tmpDF["isSS"]    = tmpDF.loc[:,"SS"] & ~(tmpDF.loc[:,"is10x"].astype('bool'))


    tmpDF["Method"] = "Unknown"
    for tech in tmpDF.columns.values[5:-1] :
        tmpDF.loc[tmpDF.loc[:,tech],"Method"]  = tech

    
    uLCP['Method'] = tmpDF.Method
    uLCP = uLCP.loc[:,["LCP","Method"]]
    df = df.merge(uLCP , on = "LCP")

    

    return(df, unknownTechs)    


def saveAsFiles(df, outpath = "MetaData",species = "mouse") : 
    # make the project directory.
    # split the dataframe by project
    by_usrp = [y for x, y in df.groupby('Project', as_index=False)]

    for i in range(len(by_usrp)) : 
        # does the file exist? 
        srp = by_usrp[i].Project.unique()[0]
        fout = outpath+"/Projects/"+species+"/"+srp+"_MetaData_TP.tsv"

        try: 
            os.makedirs(outpath+"/Projects/"+species)
        except: 
            pass

        by_usrp[i].to_csv(fout , sep="\t" , mode= "a" ,index=False, header= not os.path.exists(fout) )
    

def main( metaDirec ,species, getAll = False, getRecent=True ) : 

    # <species>_allMetaData can probably be replaced with '<species>_recentMetaData' for data pulled in the last x days. 


    # get the most recent one.
    if getAll : 
        dfpath = metaDirec +"/Projects/"+species+"_allMetaData.tsv"
        
    elif getRecent :
        dfpaths = glob.glob(metaDirec +"/Projects/"+species+"_*"  )
        dates = [path.split("_")[-1] for path in dfpaths]
        dates = [date.split(".")[0] for date in dates]
        
        # which is the most recent file?
        ind = dates.index(max(dates))
        dfpath = dfpaths[ind]
    

    # assert( os.path.exists(metaDirec +"/Projects/"+species+"_allMetaData.tsv")) 
    assert( os.path.exists(dfpath)) 

    # this may take some time to read. consider batch reads. use skiprows + nrows
    df = pd.read_csv(dfpath,sep="\t")
    df = df.loc[df.Status == 'UIDFetched' ,:]

    
    df , unknownLCP = predictTechFromLCP(df)
    df.Status = "TechPredicted"
    df.to_csv(dfpath,sep="\t",mode ="w",index=False)

    unknownLCP.to_csv(metaDirec+"/unknownLCPs.tsv",sep="\t" ,mode="a", index=False ,header=False )

    # need to save to project metadata... (individual files)


    saveAsFiles(df, outpath=metaDirec, species = species)

main(metaDirec, species, getAll=False, getRecent=True)








# if saveIt : 
#     # save the df.
#     tsvpath= h5path.replace("hdf5","tsv")
#     df.to_csv(tsvpath,sep= "\t", index=False) 

#     try: 
#         unknownLCPFile = tsvpath.replace(".tsv","_unknownTech.tsv")
#         unknownLCP.to_csv(unknownLCPFile,sep= "\t", index=False) 
#     except: 
#         pass

#     # save the method in the hdf5 file

#     with h5py.File(h5path,"a")  as f : 
#         for i in range(df.shape[0]) : 
#             ky = "/".join(df.loc[i,["Project","Experiment"] ].values)+"/"
#             try: 
#                 f[ky].create_dataset("Method", data= df.loc[i,"Method"] )
#                 print('creating ds for '+ ky)
#             except OSError as err: 
#                 # print(str(err) + " '" + ky +"'")
#                 pass
