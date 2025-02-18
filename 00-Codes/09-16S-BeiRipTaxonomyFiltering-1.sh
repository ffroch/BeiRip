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
# Filter ASV-table based on "undesirable" taxonomy
qiime taxa filter-table \
  --i-table 04-asvtable/01-raw/table.raw.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --p-exclude chloroplast,mitochondria,eukaryota \
  --o-filtered-table 04-asvtable/02-taxfilt/table.taxfilt.qza
########################################################################################################
# remove Bos taurus mitochondria and homo sapiens mitochondria
# export representative sequences
qiime tools export \
  --input-path 04-asvtable/01-raw/rep-seqs.raw.qza \
  --output-path 04-asvtable/01-raw/exported-repseqs
# blast representative sequences against human and bovine mitochondria
blastn -query 04-asvtable/01-raw/exported-repseqs/dna-sequences.fasta \
       -db /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/05-qiime/mitochondrialDB \
       -out mitochondrialhumansequences.txt \
       -outfmt "6 qseqid qlen sseqid evalue length pident"
# create a text with a column that has "FeatureID" in the first row and the first column of the txt file mitochondrialhumansequences.txt 
cut -f1 mitochondrialhumansequences.txt > featurestoremove.txt
echo "FeatureID" > featurestoremovefinal.txt
cat featurestoremove.txt >> featurestoremovefinal.txt
# filter based on featurestoremovefinal.txt
qiime feature-table filter-features \
  --i-table 04-asvtable/02-taxfilt/table.taxfilt.qza \
  --m-metadata-file featurestoremovefinal.txt \
  --p-exclude-ids \
  --o-filtered-table 04-asvtable/02-taxfilt/table.taxfilt.wobh.qza
# make barplots with filtered data
qiime taxa barplot \
  --i-table 04-asvtable/02-taxfilt/table.taxfilt.wobh.qza \
  --i-taxonomy 06-taxonomy/taxonomy.raw.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/02-taxa-bar-plots-wobh.qzv
qiime tools export \
  --input-path 07-prelimtaxabarplots/02-taxa-bar-plots.qzv \
  --output-path 07-prelimtaxabarplots/02-exported-taxa-bar-plots-wobh
########################################################################################################
conda deactivate
