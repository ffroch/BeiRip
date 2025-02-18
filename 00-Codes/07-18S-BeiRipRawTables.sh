#! /bin/bash


########################################################################################################
# Quality control of the adapter removed data
########################################################################################################
prodir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
cd $prodir/05-qiime/HTS_18S
########################################################################################################
# Initialize Conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# Activate the environment
conda activate qiime2-amplicon-2023.9 
########################################################################################################
# import data to qiime2
mkdir 04-asvtable
mkdir 04-asvtable/01-raw
# move selected truncation length files to 04-asvtable/01-raw
cp 05-trimmingoptions/table_270_210.qza 04-asvtable/01-raw/table.raw.qza
cp 05-trimmingoptions/rep-seqs_270_210.qza 04-asvtable/01-raw/rep-seqs.raw.qza
# Make visualization files
qiime feature-table summarize \
  --i-table 04-asvtable/01-raw/table.raw.qza \
  --o-visualization 04-asvtable/01-raw/table.raw.qzv \
  --m-sample-metadata-file $prodir/01-metadata/samplemeta.tsv
qiime feature-table tabulate-seqs \
  --i-data 04-asvtable/01-raw/rep-seqs.raw.qza \
  --o-visualization 04-asvtable/01-raw/rep-seqs.raw.qzv
########################################################################################################
conda deactivate
