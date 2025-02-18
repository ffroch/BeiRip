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
mkdir 09-alpha-rarefaction
# inspect 04-asvtable/table.qzv to select the optimal --p-max-depth value so you can "zoom in" the "problematic" samples
qiime diversity alpha-rarefaction \
  --i-table  04-asvtable/03-final/table.final.qza \
  --i-phylogeny 08-phylogenetictreemafft/rooted-tree.qza \
  --p-max-depth 3000 \
  --p-steps 30 \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 09-alpha-rarefaction/alpha-rarefaction.qzv \
  --p-metrics chao1 \
  --p-metrics goods_coverage \
  --p-metrics observed_features \
  --p-metrics shannon \
  --p-metrics simpson \
  --p-metrics simpson_e
########################################################################################################
conda deactivate
