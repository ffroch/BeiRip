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
# Filter taxonomy
########################################################################################################
# Step 3: filter out everything beside fungi
# Step 3a: Eukaryota without any other taxonomic information get the term "Unassigned" added
# That needed a workaraound
########################################################################################################
# features with deeper classification contain the term "d__Eukaryota;" they will be replaced by "randomplaceholder"
sed 's/d__Eukaryota;/randomplaceholder/' 06-taxonomy/02-taxonomymamfilt/taxonomy.tsv > 06-taxonomy/03-modified_taxonomy.tsv
# features classified as Eukaryota without deeper classification contain the term "d__Eukaryota" they will be replaced by "Unassigned"
sed 's/d__Eukaryota/Unassigned/' 06-taxonomy/03-modified_taxonomy.tsv > 06-taxonomy/03-modified_taxonomy1.tsv
# "randomplaceholder" will be replaced by "d__Eukaryota;" to get the original classification back
sed 's/randomplaceholder/d__Eukaryota;/' 06-taxonomy/03-modified_taxonomy1.tsv > 06-taxonomy/03-modified_taxonomy.tsv
# delete intermediate file
rm 06-taxonomy/03-modified_taxonomy1.tsv
# create a modified taxonomy qza
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat \
  --input-path 06-taxonomy/03-modified_taxonomy.tsv \
  --output-path 06-taxonomy/03-taxonomy_modified.qza
########################################################################################################
# Step3b: remove all unassigned ASVs, bacteria, archaea, apicomplexa, Foraminifera, and other undesired taxa
mkdir 04-asvtable/03-final
qiime taxa filter-table \
  --i-table 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --i-taxonomy 06-taxonomy/03-taxonomy_modified.qza \
  --p-exclude bacteria,unassigned,apicomplexa,archaea,foraminifera,embryophyta,thecofilosea,chrysophyceae,glissomonadida \
  --o-filtered-table 04-asvtable/03-final/table.clean.qza
# make preliminary taxa barplots
qiime taxa barplot \
  --i-table 04-asvtable/03-final/table.clean.qza \
  --i-taxonomy 06-taxonomy/03-taxonomy_modified.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/03-taxa-bar-plots-filtered.qzv
# export feature table
qiime feature-table summarize \
  --i-table 04-asvtable/03-final/table.clean.qza \
  --o-visualization 04-asvtable/03-final/table.clean.qzv \
  --m-sample-metadata-file $prodir/01-metadata/samplemeta.tsv
########################################################################################################
# Remove control samples based on metadata
qiime feature-table filter-samples \
  --i-table 04-asvtable/03-final/table.clean.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --p-where "[samplecontrol]='control'" \
  --p-exclude-ids --o-filtered-table 04-asvtable/03-final/table.clean.no_ntc.qza
########################################################################################################
# Remove features (ASVs) with 0 counts (or the number you want)
qiime feature-table filter-samples \
  --i-table 04-asvtable/03-final/table.clean.no_ntc.qza \
  --p-min-frequency 1 \
  --o-filtered-table 04-asvtable/03-final/table.final.qza
qiime feature-table summarize \
  --i-table 04-asvtable/03-final/table.final.qza \
  --o-visualization 04-asvtable/03-final/table-final.qzv \
  --m-sample-metadata-file $prodir/01-metadata/samplemeta.tsv
########################################################################################################
# Remove reads that are not in the ASV table anymore (not mandatory, but will save time and resources for downstream analyses)
qiime feature-table filter-seqs \
  --i-data 04-asvtable/01-raw/rep-seqs.raw.qza \
  --i-table 04-asvtable/03-final/table.final.qza \
  --o-filtered-data 04-asvtable/03-final/rep-seqs.final.qza
# export rep-seqs
qiime tools export \
  --input-path 04-asvtable/03-final/rep-seqs.final.qza \
  --output-path 04-asvtable/03-final/exportedrepseqs.final
########################################################################################################
# check final table
qiime taxa barplot \
  --i-table 04-asvtable/03-final/table.final.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/04-taxa-bar-plots-final.qzv
qiime tools export \
  --input-path 07-prelimtaxabarplots/04-taxa-bar-plots-final.qzv \
  --output-path 07-prelimtaxabarplots/04-exported-taxa-bar-plots-final
########################################################################################################
conda deactivate
