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
# Export ASV table into a txt/tsv/csv file that can be read in R
qiime tools export \
  --input-path 04-asvtable/02-taxfilt/table.taxfilt.wobh.qza \
  --output-path 04-asvtable/02-taxfilt/exported-table-taxfilt-wobh
biom convert \
  -i 04-asvtable/02-taxfilt/exported-table-taxfilt-wobh/feature-table.biom \
  -o 04-asvtable/02-taxfilt/exported-table-taxfilt-wobh/feature-table.txt \
  --to-tsv
sed -i '1d' 04-asvtable/02-taxfilt/exported-table-taxfilt-wobh/feature-table.txt
sed -i "s/^#OTU ID/OTU_ID/" 04-asvtable/02-taxfilt/exported-table-taxfilt-wobh/feature-table.txt
########################################################################################################
conda deactivate
