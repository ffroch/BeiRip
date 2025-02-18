#! /bin/bash
########################################################################################################
# Adapter removal
########################################################################################################
prodir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
########################################################################################################
# Initialize Conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# Activate the environment
conda activate trim
########################################################################################################
# 16S
########################################################################################################
cd $prodir/05-qiime/HTS_16S
# Adapters removal using trimmomatic
for i in $(ls -1 $prodir/02-rawdata/HTS_16S/fastq/*_R1_* | xargs -n 1 basename | sed "s/_R1_001.fastq.gz//"); do
  trimmomatic PE \
    -threads 24 \
    -phred33 $prodir/02-rawdata/HTS_16S/fastq/${i}_R1_001.fastq.gz \
      $prodir/02-rawdata/HTS_16S/fastq/${i}_R2_001.fastq.gz \
      02-noadapt/${i}_R1_001.fastq.gz /dev/null 02-noadapt/${i}_R2_001.fastq.gz \
      /dev/null ILLUMINACLIP:/data/Unit_LMM/databases/seq-adapters.fasta:1:30:11
done
########################################################################################################
# 18S
########################################################################################################
cd $prodir/05-qiime/HTS_18S
# Adapters removal using trimmomatic
conda activate trim
for i in $(ls -1 $prodir/02-rawdata/HTS_18S/fastq/*_R1_* | xargs -n 1 basename | sed "s/_R1_001.fastq.gz//"); do
  trimmomatic PE \
    -threads 24 \
    -phred33 $prodir/02-rawdata/HTS_18S/fastq/${i}_R1_001.fastq.gz \
      $prodir/02-rawdata/HTS_18S/fastq/${i}_R2_001.fastq.gz \
      02-noadapt/${i}_R1_001.fastq.gz /dev/null 02-noadapt/${i}_R2_001.fastq.gz \
      /dev/null ILLUMINACLIP:/data/Unit_LMM/databases/seq-adapters.fasta:1:30:11
done
########################################################################################################
conda deactivate
########################################################################################################