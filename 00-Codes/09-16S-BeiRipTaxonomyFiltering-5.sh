#! /bin/bash


########################################################################################################
# Re-import data to qiime after decontam
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
# Import the output from decontam into Q2
## Make sure that the header starts with "#OTU ID", for instance:
sed "1s/^/#OTU ID/" 04-asvtable/decontam/feature-table.decontam.txt

biom convert \
  -i 04-asvtable/decontam/feature-table.decontam.txt \
  -o 04-asvtable/decontam/feature-table.decontam.biom \
  --to-hdf5
qiime tools import \
  --input-path 04-asvtable/decontam/feature-table.decontam.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 04-asvtable/03-decontam/table.taxfilt.decontam.qza

# Remove control samples based on metadata
qiime feature-table filter-samples \
  --i-table 04-asvtable/03-decontam/table.taxfilt.decontam.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --p-where "[samplecontrol]='control'" \
  --p-exclude-ids \
  --o-filtered-table 04-asvtable/03-decontam/table.taxfilt.decontam.no_ntc.qza

# Remove features (ASVs) with 0 counts (or the number you want)
qiime feature-table filter-samples \
  --i-table 04-asvtable/03-decontam/table.taxfilt.decontam.no_ntc.qza \
  --p-min-frequency 1 \
  --o-filtered-table 04-asvtable/04-final/table.qza
qiime feature-table summarize \
  --i-table 04-asvtable/04-final/table.qza \
  --o-visualization 04-asvtable/04-final/table.qzv \\
  --m-sample-metadata-file $prodir/01-metadata/samplemeta.tsv

# Remove reads that are not in the ASV table anymore (not mandatory, but will save time and resources for downstream analyses)
qiime feature-table filter-seqs \
  --i-data 04-asvtable/01-raw/rep-seqs.raw.qza \
  --i-table 04-asvtable/04-final/table.qza \
  --o-filtered-data 04-asvtable/04-final/rep-seqs.qza

# check final table
qiime taxa barplot \
  --i-table 04-asvtable/04-final/table.qza \
  --i-taxonomy 06-taxonomy/01-taxonomy.raw.qza \
  --m-metadata-file $prodir/01-metadata/samplemeta.tsv \
  --o-visualization 07-prelimtaxabarplots/03-taxa-bar-plots-postdecontam.qzv
qiime tools export \
  --input-path 07-prelimtaxabarplots/taxa-bar-plots-postdecontam.qzv \
  --output-path 07-prelimtaxabarplots/03-exported-taxa-bar-plots-postdecontam
########################################################################################################
# make blast db
mkdir -p 08-blastdb
# export repseqs from the final ASV-table
qiime tools export \
  --input-path 04-asvtable/04-final/rep-seqs.qza \
  --output-path 08-blastdb
cd 08-blastdb
# make blast db
makeblastdb \
  -in dna-sequences.fasta \
  -title rep-seqs \
  -dbtype nucl \
  -hash_index
########################################################################################################
# blast isolates against sequencing data
blastn \
  -query /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/03-Routputfiles/C1_seqstrimclean_35_15_100bp_bac.fa \
  -db dna-sequences.fasta \
  -out comparison.blastout \
  -num_alignments 25 \
  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads 100
cd ../..
########################################################################################################
conda deactivate
