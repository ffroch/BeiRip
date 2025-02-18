#! /bin/bash


########################################################################################################
# Quality control of the adapter removed data
########################################################################################################
prodir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
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
mkdir -p 01-qc/fastqc-no-adapt 01-qc/multiqc-no-adapt
# Quality Control of the adapter removed data
for i in $(ls -1 02-noadapt/*fastq.gz | xargs -n 1 basename); do
  fastqc 02-noadapt/${i} -o 01-qc/fastqc-no-adapt -t 24
done
# summarize the results of the fastqc analysis
multiqc 01-qc/fastqc-no-adapt/*_R1* -o 01-qc/multiqc-no-adapt -i R1
multiqc 01-qc/fastqc-no-adapt/*_R2* -o 01-qc/multiqc-no-adapt -i R2
########################################################################################################
# 18S
########################################################################################################
cd $prodir/05-qiime/HTS_18S
mkdir -p 01-qc/fastqc-no-adapt 01-qc/multiqc-no-adapt
# Quality Control of the adapter removed data
for i in $(ls -1 02-noadapt/*fastq.gz | xargs -n 1 basename); do
  fastqc 02-noadapt/${i} -o 01-qc/fastqc-no-adapt -t 24
done
# summarize the results of the fastqc analysis
multiqc 01-qc/fastqc-no-adapt/*_R1* -o 01-qc/multiqc-no-adapt -i R1
multiqc 01-qc/fastqc-no-adapt/*_R2* -o 01-qc/multiqc-no-adapt -i R2
########################################################################################################
conda deactivate
########################################################################################################