#! /bin/bash


########################################################################################################
# Quality control of the adapter removed data
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
# import data to qiime2
mkdir 03-qiimeraw
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 02-noadapt \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path 03-qiimeraw/PE-reads.qza
# inspect meta data file
qiime tools inspect-metadata $prodir/01-metadata/samplemeta.tsv #IMPORTANT save tsv as UTF-8!!!
# make a quality plot of the reads with QIIME2
qiime demux summarize \
  --i-data 03-qiimeraw/PE-reads.qza \
  --o-visualization 03-qiimeraw/PE-reads.qzv
########################################################################################################
# run filtering with different truncation length combinations to evaluate are the best for your data
# I used trim-left-f 18 and trim-left-r 18, because that is the length of the primers used
mkdir 05-trimmingoptions
for i in 250 260 270 280 290; do
for j in 210 220 230 240 250 260; do
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 03-qiimeraw/PE-reads.qza \
  --p-trim-left-f 18 \
  --p-trim-left-r 18 \
  --p-trunc-len-f ${i} \
  --p-trunc-len-r ${j} \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-pooling-method 'pseudo' \
  --p-chimera-method 'consensus' \
  --o-representative-sequences 05-trimmingoptions/rep-seqs_${i}_${j}.qza \
  --o-table 05-trimmingoptions/table_${i}_${j}.qza \
  --o-denoising-stats 05-trimmingoptions/stats-dada2_${i}_${j}.qza \
  --p-n-threads 100
qiime metadata tabulate \
  --m-input-file 05-trimmingoptions/stats-dada2_${i}_${j}.qza \
  --o-visualization 05-trimmingoptions/stats-dada2_${i}_${j}.qzv
qiime tools export \
  --input-path 05-trimmingoptions/stats-dada2_${i}_${j}.qzv \
  --output-path 05-trimmingoptions/stats-dada2_${i}_${j}
done
done
########################################################################################################
conda deactivate
