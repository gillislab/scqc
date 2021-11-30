#!/usr/bin/env bash
set -e
# install conda
if ! command -v conda&> /dev/null; then
	echo "installing miniconda..."
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
	bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b
	rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh
	conda init 
	echo "miniconda installed. restart terminal."
	exit 0
else
	echo "miniconda installed already."
fi
# be sure to restart terminal to allow conda to start

# to export environment to yml
# conda env export --name scqc > scqc.yml
# to create conda from yml file. Otherwise, skip and install manually
# conda env create --file scqc.yml

conda create -y -n scqc python=3.8
conda activate scqc

# make sure we have conda-forge
conda config --add channels conda-forge

# update conda if necessary
conda update -y -n base -c defaults conda

conda install -y  pandas ipython requests scipy wget bottleneck tbb
conda install -y -c bioconda star
conda install -y -c conda-forge leidenalg scanpy scikit-misc

#download sratoolkit, link binaries within conda environment
cd  $CONDA_PREFIX

# linux/centos
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-centos_linux64.tar.gz
tar -xvzf sratoolkit.2.11.0-centos_linux64.tar.gz
cd bin ; ln -s ../sratoolkit.2.11.0-centos_linux64/bin/* ./

# macos
# wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-mac64.tar.gz
# tar -xvzf sratoolkit.2.11.0-mac64.tar.gz
# cd bin
#ln -s ../sratoolkit.2.11.0-mac64/bin/* ./

#  R no longer needed....
#
# install R and dependencies on conda
#conda config --set channel_priority strict
# required for R 4.0
#conda install -c conda-forge r-base 
#conda install -c conda-forge r-essentials
#conda install -c conda-forge r-cairo 
#conda install -c conda-forge r-ragg


# R pacakges to install
#####  R no longer needed with Python metamarkers...
#R -e "install.packages('devtools',repos='http://cran.us.r-project.org')"
#R -e "install.packages('BiocManager',repos='http://cran.us.r-project.org')"
#R -e "install.packages('optparse',repos='http://cran.us.r-project.org')"
#R -e "install.packages('Matrix',repos='http://cran.us.r-project.org')"
#R -e "remove.packages('textshaping')"
#R -e "BiocManager::install('scRNAseq')"
#R -e "BiocManager::install('rhdf5')"
#R -e "devtools::install_github('gillislab/MetaMarkers')"
