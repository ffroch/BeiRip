#! /bin/bash


########################################################################################################
# Make Phalogenetic Tree
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
mkdir 08-phylogenetictreemafft
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 04-asvtable/03-final/rep-seqs.final.qza \
  --o-alignment 08-phylogenetictreemafft/aligned-rep-seqs.qza \
  --o-masked-alignment 08-phylogenetictreemafft/masked-aligned-rep-seqs.qza \
  --o-tree 08-phylogenetictreemafft/unrooted-tree.qza \
  --o-rooted-tree 08-phylogenetictreemafft/rooted-tree.qza \
  --p-n-threads 100
########################################################################################################
conda deactivate
