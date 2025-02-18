#! /bin/bash


########################################################################################################
# Quality control of the adapter removed data
########################################################################################################
prodir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
cd $prodir/05-qiime/HTS_16S
########################################################################################################
# Initialize Conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
# Activate the environment
conda activate qiime2-amplicon-2023.9 
########################################################################################################
# import data to qiime2
mkdir 06-taxonomy
# taxonomic classification
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Unit_LMM/selberherr-group/roch/00_databases/full-length-uniform-classifier.qza \
  --i-reads 04-asvtable/01-raw/rep-seqs.raw.qza \
  --o-classification 06-taxonomy/01-taxonomy.raw.qza \
  --p-n-jobs 100
# export taxonomy
qiime tools export \
  --input-path 06-taxonomy/01-taxonomy.raw.qza \
  --output-path 06-taxonomy/01-exported-taxonomy-raw
# make preliminary taxa barplots
mkdir 07-prelimtaxabarplots
qiime taxa barplot \
  --i-table 04-asvtable/01-raw/table.raw.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/01-taxa-bar-plots.qzv
qiime tools export \
  --input-path 07-prelimtaxabarplots/01-taxa-bar-plots.qzv \
  --output-path 07-prelimtaxabarplots/01-exported-taxa-bar-plots
########################################################################################################
conda deactivate
