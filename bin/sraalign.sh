#!/bin/bash
# Let's first re-generate the genome using an up-to-date STAR
# ./STAR  --runThreadN 24 --runMode genomeGenerate --genomeDir /data/johlee/forHamsini/yeastGenome --genomeFastaFiles /data/johlee/forHamsini/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa --sjdbGTFfile /data/johlee/forHamsini/Saccharomyces_cerevisiae.R64-1-1.93.gtf --sjdbOverhang 100
# only needs to be done once.
dirname="/home/johlee/QC_project"
cd $dirname
genDir="/data/genomes/GRCm38"
SRA="SRP185852"   
# get the runs associated with this SRA ID
wget -O ${SRA}'_metaData.csv' "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="${SRA}
# list of SRA accessions that need to be realigned i.e. SRP##### 
cut -f1 -d, ${SRA}"_metaData.csv" > ${SRA}"_runs.csv"
sed -i 1d ${SRA}"_runs.csv"         # removes colname
nCoreSTAR=20  

getFastq() {
    local Run=$1  # first position
    link="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    first6="$(cut -b 1-6 <<<"$Run")"
    last="${Run: -1}"
    L=$(echo -n $Run | wc -m)
    if [ $L = 10 ]  ;then # 7 digits
        link="${link}${first6}/00${last}/${Run}/*"
        wget -nv -nc  "$link"
    elif [ $L = 9 ]; then # 6 digits
        link="${link}${first6}/${Run}/*"
        wget -nv -nc  "$link"
    fi 
}


runSTAR() {
    local Run=$1 
    local nCore=$2 
    local genDir=$3 
    local filein=$4
    ./STAR --runThreadN $nCore --genomeDir $genDir --readFilesIn $filein --readFilesCommand zcat --outSAMtype None --outFileNamePrefix $Run --quantMode GeneCounts
    if test -f "${Run}ReadsPerGene.out.tab" ; then
        rm ${Run}*.fastq.gz
    fi
}

SRAfile=${SRA}"_runs.csv"  # runs within the SRA    
while IFS= read -r Run
do
    echo $Run
    # get the fastq file(s)
    getFastq $Run 
    # is it paired?
    filein=$Run".fastq.gz"
    filesin_paired=$Run"_1.fastq.gz"
    if test -f $filein; then
        filein=$filein
    elif test -f $filesin_paired; then
        filein="${filesin_paired} ${Run}_2.fastq.gz"
    fi
    # runs STAR using default parameters
    runSTAR $Run $nCoreSTAR $genDir $filein 
done < $SRAfile

