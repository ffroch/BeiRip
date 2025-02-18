#! /bin/bash

########################################################################################################
# Quality control of the raw data
########################################################################################################
prodir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
# go to working directory
cd $prodir
# create directory called 05-qiime
mkdir 05-qiime
cd 05-qiime
mkdir HTS_16S HTS_18S
########################################################################################################
# Initialize Conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# Activate the environment
conda activate QC
########################################################################################################
# 16S
########################################################################################################
cd $prodir/05-qiime/HTS_16S
# create directory called 01-qc	and 02-noadapt
mkdir 01-qc
# create subdirectories for fastqc and multiqc
mkdir -p 01-qc/fastqc-raw 01-qc/multiqc-raw

# Quality Control of the raw data
for i in $(ls -1 $prodir/02-rawdata/HTS_16S/fastq/*fastq.gz | xargs -n 1 basename); do
  fastqc $prodir/02-rawdata/HTS_16S/fastq/${i} -o 01-qc/fastqc-raw -t 24
done
multiqc 01-qc/fastqc-raw/*_R1* -o 01-qc/multiqc-raw -i R1
multiqc 01-qc/fastqc-raw/*_R2* -o 01-qc/multiqc-raw -i R2

# renaming files to match import format (Casava 1.8 paired-end demultiplexed fastq)
cd $prodir/02-rawdata/HTS_16S/fastq
# remove _LBFVB.filt from all the file names in this directory
for file in *"_LBFVB.filt"*; do
    # Use parameter expansion to remove the term
    newname="${file/_LBFVB.filt/}"
    mv "$file" "$newname"
done
# remove _L8DG4.filt from all the file names in this directory
for file in *"_L8DG4.filt"*; do
    # Use parameter expansion to remove the term
    newname="${file/_L8DG4.filt/}"
    mv "$file" "$newname"
done
########################################################################################################
# 18S
########################################################################################################
cd $prodir/05-qiime/HTS_18S
# create directory called 01-qc	and 02-noadapt
mkdir 01-qc 
# create subdirectories for fastqc and multiqc
mkdir -p 01-qc/fastqc-raw 01-qc/multiqc-raw

# Quality Control of the raw data
for i in $(ls -1 $prodir/02-rawdata/HTS_18S/fastq/*fastq.gz | xargs -n 1 basename); do
  fastqc $prodir/02-rawdata/HTS_18S/fastq/${i} -o 01-qc/fastqc-raw -t 24
done
multiqc 01-qc/fastqc-raw/*_R1* -o 01-qc/multiqc-raw -i R1
multiqc 01-qc/fastqc-raw/*_R2* -o 01-qc/multiqc-raw -i R2

# renaming files to match import format (Casava 1.8 paired-end demultiplexed fastq)
cd $prodir/02-rawdata/HTS_18S/fastq
# remove _LBFVB.filt from all the file names in this directory
for file in *"_LBFVB.filt"*; do
    # Use parameter expansion to remove the term
    newname="${file/_LBFVB.filt/}"
    mv "$file" "$newname"
done
# remove _LCP57.filt from all the file names in this directory
for file in *"_LCP57.filt"*; do
    # Use parameter expansion to remove the term
    newname="${file/_LCP57.filt/}"
    mv "$file" "$newname"
done
########################################################################################################
# Deactivate the environment
conda deactivate