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
qiime tools export \
  --input-path 06-taxonomy/01-taxonomy.raw.qza \
  --output-path 06-taxonomy/01-exported-taxonomy-raw
# export feature table
qiime tools export \
  --input-path 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --output-path 04-asvtable/02-mamfilt/exported-table-mamfilt
# convert biom to tsv
biom convert \
  -i 04-asvtable/02-mamfilt/exported-table-mamfilt/feature-table.biom \
  -o 04-asvtable/02-mamfilt/exported-table-mamfilt/feature-table.tsv \
  --to-tsv
# remove taxonomy based on filtered ASVs
grep -Ff <(cut -f1 04-asvtable/02-mamfilt/exported-table-mamfilt/feature-table.tsv) 06-taxonomy/01-exported-taxonomy-raw/taxonomy.tsv > 06-taxonomy/02-filtered_taxonomy.tsv
# make taxonomy qza
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 06-taxonomy/02-filtered_taxonomy.tsv \
  --output-path 06-taxonomy/02-taxonomy.mamfilt.qza
########################################################################################################
# export relevant files
#######################################################################################################
# filter rep-seqs from the lates ASV-table
qiime feature-table filter-seqs \
  --i-data 04-asvtable/01-raw/rep-seqs.raw.qza \
  --i-table 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --o-filtered-data 04-asvtable/02-mamfilt/rep-seqs.mamfilt.qza
# export feature table 
qiime tools export \
  --input-path 04-asvtable/02-mamfilt/table.mamfilt.qza \
  --output-path 08-blastdb/table.mamfilt
# convert biom to tsv
biom convert \
  -i 08-blastdb/table.mamfilt/feature-table.biom \
  -o 08-blastdb/table.mamfilt/feature-table.txt \
  --to-tsv
# remove first line
tail -n+2 08-blastdb/table.mamfilt/feature-table.txt | wc -l
# export rep-seqs
qiime tools export \
  --input-path 04-asvtable/02-mamfilt/rep-seqs.mamfilt.qza \
  --output-path 04-asvtable/02-mamfilt/repseqsmamfilt
# export taxonomy
qiime tools export \
  --input-path 06-taxonomy/02-taxonomy.mamfilt.qza  \
  --output-path 06-taxonomy/02-taxonomymamfilt
# copy folder repseqsmamfilt to blastdb
cp -r 04-asvtable/02-mamfilt/repseqsmamfilt 08-blastdb/
########################################################################################################
# make blast db
########################################################################################################
cd 08-blastdb/repseqsmamfilt/
# make blast db
makeblastdb \
  -in dna-sequences.fasta \
  -title rep-seqs \
  -dbtype nucl \
  -hash_index
########################################################################################################
# blast isolates against sequencing data
########################################################################################################
blastn \
  -query /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/03-Routputfiles/C2_seqstrimclean_35_15_100bp_fun.fa \
  -db dna-sequences.fasta \
  -out ../comparison.mamfilt.blastout \
  -num_alignments 25 \
  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads 100
########################################################################################################
cd ../..
conda deactivate
