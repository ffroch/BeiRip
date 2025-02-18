#! /bin/bash


########################################################################################################
# Prepare data for R
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
# export ASV table
qiime tools export \
  --input-path 04-asvtable/04-final/table.qza \
  --output-path 04-asvtable/04-final/exported-table
biom convert \
  -i 04-asvtable/04-final/exported-table/feature-table.biom \
  -o 04-asvtable/04-final/exported-table/feature-table.txt \
  --to-tsv
sed -i '1d' 04-asvtable/04-final/exported-table/feature-table.txt
sed -i "s/^#OTU ID/OTU_ID/" 04-asvtable/04-final/exported-table/feature-table.txt
# export phylogenetic tree
qiime tools export \
  --input-path 08-phylogenetictreemafft/rooted-tree.qza \
  --output-path 08-phylogenetictreemafft/exported-rooted-tree
########################################################################################################
conda deactivate
