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
# Step 1: filter out all ASVs that blast to Bos taurus
# I added that steps, since many of the ASVs were in general unassigned or unassigned Eukaryota after classification, but blast showed that they were associated with Bos taurus or other Artiodactyla. With this step I hope to reduce those ASVs a bit. 
# export rep-seqs
qiime tools export \
  --input-path 04-asvtable/01-raw/rep-seqs.raw.qza \
  --output-path 04-asvtable/01-raw/exported-repseqs
# blast rep-seqs against Bos taurus genome
blastn \
  -query 04-asvtable/01-raw/exported-repseqs/dna-sequences.fasta \
  -db /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/05-qiime/bos_taurus_DB \
  -out bostaurushits.txt \
  -outfmt "6 qseqid qlen sseqid evalue length pident" \
  -num_threads 100
# keep only hit with > 240 bp length
awk '$5 >= 240 { print }' bostaurushits.txt > bostaurushits240.txt
# keep only 100% pident hits
awk '$6 >= 97.00 { print }' bostaurushits240.txt > bostaurushits97.txt
# remove duplicates
sort -k1,1 -u bostaurushits97.txt > bostaurushits97uniq.txt
# i check this file manually and compared it with the taxonomy file. Nothing that was classified as fungus was listed.
# create a text with a column that has "FeatureID" in the first row and the first column of the txt file mitochondrialhumansequences.txt 
cut -f1 bostaurushits97uniq.txt > bostaurushitstoremove.txt
echo "FeatureID" > bostaurushitstoremovefinal.txt
cat bostaurushitstoremove.txt >> bostaurushitstoremovefinal.txt
# filter based on featurestoremovefinal.txt
qiime feature-table filter-features \
 --i-table 04-asvtable/01-raw/table.raw.qza \
 --m-metadata-file bostaurushitstoremovefinal.txt \
 --p-exclude-ids \
 --o-filtered-table 04-asvtable/02-mamfilt/01-table.wobt.qza
########################################################################################################
conda deactivate
