#! /bin/bash
########################################################################################################
# FILE PREPARATION AND ORGANIZATION    #################################################################
########################################################################################################
# this project had to kind of sequencing data. 16S and 18S so the data was organized in two folders

# go to project directory
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/02-rawdata
# make two directories for 16S and 18S analysis
mkdir HTS_16S HTS_18S

# unzip file from Baseclear
unzip 158031.zip

# the following numbers are extracted from the Basclear Reports called 158031_Next-Generation-Sequencing-project_report.pdf and 158031_Next-Generation-Sequencing-project_report_v2.pdf, both saved in the folder 02-rawdata
# move all files containing 109022, 109023, 109699, 110021, 110097, and 110098 from the subdirectory raw-sequences to the subdirectory HTS_16S
for i in 109022 109023 109699 110021 110097 110098; do
  mv raw_sequences/*${i}* HTS_16S/
done
# move all files containing 109919, 109920, and 110099 from the subdirectory raw-sequences to the subdirectory HTS_18S
for i in 109919 109920 110099 111828 111829; do
  mv raw_sequences/*${i}* HTS_18S/
done
# remove the subdirectory raw-sequences
rm -r raw_sequences
########################################################################################################
# move company generated fastqc files and fastq files into separate folders for 16S and 18S analysis  
cd HTS_16S
mkdir companyfastqc fastq
# move all files containing fastqc into the folder companyfastqc
for i in $(ls -1 *fastqc*); do
  mv $i companyfastqc/
done
# move all files containing fastq into the folder fastq
for i in $(ls -1 *fastq*); do
  mv $i fastq/
done
cd ../HTS_18S
mkdir companyfastqc fastq
# move all files containing fastqc into the folder companyfastqc
for i in $(ls -1 *fastqc*); do
  mv $i companyfastqc/
done
# move all files containing fastq into the folder fastq
for i in $(ls -1 *fastq*); do
  mv $i fastq/
done
########################################################################################################