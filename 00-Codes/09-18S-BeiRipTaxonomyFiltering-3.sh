#! /bin/bash


########################################################################################################
# Filter ASV-table based on "undesirable" taxonomy
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
## Step 2: filter out all ASVs that are mammalia
# first I want to keep all non mammalian ASVs in the table. I found it interesting to see, in how many samples there were sarcocystis ASVs. Further I want to check if we have a lot of Archea in it. The bacteria I keep because I want to check if I could use them as internal standard in combination with the 16S data. 
qiime taxa filter-table \
  --i-table 04-asvtable/02-mamfilt/01-table.wobt.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --p-exclude mammalia \
  --o-filtered-table 04-asvtable/02-mamfilt/table.mamfilt.qza
# make preliminary taxa barplots
qiime taxa barplot \
  --i-table 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/02-taxa-bar-plots-mammaliafiltered.qzv
# export feature table
qiime feature-table summarize \
  --i-table 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --o-visualization 04-asvtable/02-mamfilt/table.mamfilt.qzv \
  --m-sample-metadata-file $prodir/01-metadata/samplemeta.tsv
########################################################################################################
conda deactivate
