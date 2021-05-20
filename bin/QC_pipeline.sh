#!/bin/bash


# Requires entrez-direct - No longer used
# conda install -c bioconda entrez-direct

# Requires STAR
# wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
# tar -xzf 2.7.8a.tar.gz
# cd STAR-2.7.8a/source
# make STAR
# export PATH=$PATH:$PWD/STAR-2.7.8a/bin/Linux_x86_64

# Requires sratoolkit
# wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
# tar -vxzf sratoolkit.tar.gz
# export PATH=$PATH:$PWD/sratoolkit<__extension__>/bin




# to add more species, links to the fa/gtf files must be provided in "getGenomeData"
# then use STAR to generate the genome for the new species

# no longer useing Svensson's dataset
# use: getSvenssonData $suppDirec
getSvenssonData(){
    local data_direc=$1
    data_direc=$data_direc"/svensson"   

    mkdir -p $data_direc
    
    current_date=$(date "+%Y%m%d")
    wget -nv -O $data_direc"/svensson"$current_date".tsv" http://www.nxn.se/single-cell-studies/data.tsv 
    cp $data_direc"/svensson"$current_date".tsv" $data_direc"/svenssonCurrent.tsv"
}

# only needs to be done once
# use: getWhiteLists $suppDirec
getWhiteLists(){
    local data_direc=$1
    mkdir -p $data_direc"/whiteLists"

    wget -nv -nc -O $data_direc"/whiteLists/3M-february-2018.txt.gz" https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
    yes y | gunzip $data_direc"/whiteLists/3M-february-2018.txt.gz"

    wget -nv -nc -O $data_direc"/whiteLists/737K-august-2016.txt" https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
    
    wget -nv -nc -O $data_direc"/whiteLists/737K-april-2014_rc.txt" https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt

}

# download species gtf/fasta files
# needs to be updated when adding a new species
# use: $ getGenomeData $species $data_direc
getGenomeData(){
    local species=$1
    local data_direc=$2 
    # need to add an update option
    
    echo "... ... Getting fasta and gtf files for "$species
    # links need to be specified for additional species
    if [ $species == "human" ] ; then 
        fastaLink="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz"
        gtfLink="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz"
    elif [ $species == "mouse" ] ; then 
        fastaLink="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.primary_assembly.genome.fa.gz"
        gtfLink="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz"
    else 
        echo "... ... No links provided for "$species" yet. Add them to function 'getGenomeData' "
        return
    fi 


    mkdir -p $data_direc"/genomes/"$species

    # is the fasta file already there? dont over write
    if ! test -f $data_direc"/genomes/"$species"/genome.fa" ; then 
        wget -nv -nc -O $data_direc"/genomes/"$species"/genome.fa.gz" $fastaLink
        gunzip  $data_direc"/genomes/"$species"/genome.fa.gz" 
    fi 

    if ! test -f $data_direc"/genomes/"$species"/annotation.gtf" ; then
            wget -nv -nc -O $data_direc"/genomes/"$species"/annotation.gtf.gz" $gtfLink
            gunzip -f $data_direc"/genomes/"$species"/annotation.gtf.gz"
    fi

}

# used only once per species.
# generateGenome $species $data_direc $nCore "FALSE"
generateGenome(){
    local species=$1    
    local data_direc=$2
    local nCore=$3
    local update=$4 # should this be updated?
    
    echo "... Generating the genome for"$species 
    # does the directory already exist?? shouldn't regenerate it.
    if test -d $data_direc"/genomes/"$species"/STARindices" ; then
        if [ $update == "TRUE" ]; then
            rm -r $data_direc"/genomes/"$species"/STARindices"
        else
            echo "... ... Genome already generated for "$species" using STAR. Nothing to do."  
            return 
        fi 
    fi 


    # get the data if not already available
    getGenomeData $species $data_direc
    wait

    # are the fasta files there??
    if  ! test -f $data_direc"/genomes/"$species"/annotation.gtf" ; then 
        echo "... Can't find GTF file for "$species". Exiting."
        return  
    elif ! test -f $data_direc"/genomes/"$species"/genome.fa" ; then
        echo "... Can't find FASTA file for "$species". Exiting."
        return  
    fi

    # use STAR to generate the index
    if test -f "STAR" ; then
        mkdir -p $data_direc"/genomes/"$species"/STARindices"

        echo "... Running STAR"

        path2fasta=$data_direc"/"
        ./STAR \
            --runMode genomeGenerate \
            --runThreadN $nCore \
            --genomeDir $data_direc"/genomes/"$species"/STARindices" \
            --genomeFastaFiles $data_direc"/genomes/"$species"/genome.fa" \
            --sjdbGTFfile $data_direc"/genomes/"$species"/annotation.gtf"
    fi


}



###################### Temporary config file ######################
workingDirec="/data/johlee/QC_project"
cd $workingDirec

suppDirec="SupplementData"
metaDataDirec="MetaData"
fastqDirec="FASTQ"

nCore="20" 
species="mouse"     
GSEOnly="TRUE"      # should be specified (TRUE/FALSE), default is FALSE
Species="Mouse"     # Mouse, Human, NA
Tissue="Brain"      # Brain, Blood, etc, NA
Tech="Chromium"     # Chromium, Smart-seq, SMARTer, etc, NA

mkdir -p $suppDirec
mkdir -p $metaDataDirec
mkdir -p $fastqDirec

################################################################

################ get supp data if not available ################

# get data from svensson's database
getSvenssonData $suppDirec
# get whitelists for STARsolo
getWhiteLists $suppDirec
# generate the genomes for STAR for the specified species 
generateGenome $species $suppDirec $nCore "FALSE" # If True - will update the STARindex. NOT recommended

#################################################################

####################### get svensson data #######################
Rscript parseDatasets.R  $GSEOnly $Species $Tissue $Tech $suppDirec"/svensson"
GSEList=$suppDirec'/svensson/svensson_'$Species'_'$Tissue"_"$Tech"_GSEs.tsv"


# need to verify if these are already done. 
# convert these GSE ids to SRA and get the corresponding metadata.
Rscript getMetaDataFromGSE.R $GSEList $suppDirec $metaDataDirec $nCore

# run STAR solo given the options.
# begin for loop around SRA...
SRA="SRA578634"



# downloads fastq files in parallel using sratoolkit's prefetch + fasterq-dump
filepath=$metaDataDirec"/"$SRA"/"$SRA"_metaData.tsv"
Rscript getFastqFiles.R $filepath $fastqDirec $nCore




# Builds a manifest file for each technology used in the SRA
Rscript getSTARparams.R "$fastqDirec$SRA"

# runs STAR - looks in fastqDirec for the manifest file
Rscript STARwrapper.R $fastqDirec $suppDirec $species

# samtools fastq -@ 8 SRR5799988.bam \
#     -1 SAMPLE_R1.fastq.gz \
#     -2 SAMPLE_R2.fastq.gz \
#     -0 /dev/null -s /dev/null -n