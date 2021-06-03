import pandas as pd
# import h5py

solooutDir ='/home/johlee/SCQC/STARout/SRR11604218Solo.out'

def gather_stats_from_STAR(solooutDir): 

    with open(f"{solooutDir}/Barcodes.stats" ,'r') as f :
        lines = f.readlines()
        lines =  [ " ".join(line.split() )  for line in lines ]
        lines =  [ line.split()  for line in lines ]
        barcode_stats = pd.DataFrame(lines, columns = ["Stat", "Value"])

    
    with open(f"{solooutDir}/Gene/Features.stats" ,'r') as f :
        lines = f.readlines()
        lines =  [ " ".join(line.split() )  for line in lines ]
        lines =  [ line.split()  for line in lines ]
        feature_stats = pd.DataFrame(lines, columns = ["Stat", "Value"])

    summary_stats = pd.read_csv(f"{solooutDir}/Gene/Summary.csv" ,sep=",",header=None)
    summary_stats.columns = ["Stat", "Value"]

    acc = solooutDir.split("/")[-1].split("Solo.out")[0]
    barcode_stats["Accession"] = acc
    feature_stats["Accession"] = acc
    summary_stats["Accession"] = acc

    return(barcode_stats, feature_stats , summary_stats )




