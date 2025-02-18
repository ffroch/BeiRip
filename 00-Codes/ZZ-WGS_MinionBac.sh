##################################################################################
# USED WORKFLOW FOR MINION BACTERIAL GENOME ASSEMBLY AND ANNOTATION
##################################################################################
# 1. Quality check of the raw data
conda activate Minion-2024

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/06-Minion
cd $workdir
mkdir -p 02-qc

# check the quality of the raw data
NanoPlot --summary ./01-rawdata/Minion20241024/sequencing_summary_FBA76834_901686b0_fea538ea.txt \
    -o ./02-qc/Minion20241024

NanoPlot --summary ./01-rawdata/Minion20241024/sequencing_summary_FBA76834_901686b0_fea538ea.txt \
    -o ./02-qc/Minion20241024-filtered --minlength 1000 --minq 7

##################################################################################
# 2. Convert bam files to fastq files
mkdir -p 05-fastq
for i in {01..24}; do
    samtools fastq 01-demux/*barcode$i.bam > 05-fastq/barcode$i.fastq
done

##################################################################################
# 3. Quality check of the fastq files
for i in {01..24}; do
    echo "Stats for barcode$i.fastq" >> 02-qc/all_barcodes.stats.txt
    seqkit stats 05-fastq/barcode$i.fastq >> 02-qc/all_barcodes.stats.txt
    echo -e "\n" >> 02-qc/all_barcodes.stats.txt # Add a newline for readability
done

# also check quality with NanoPlot
for i in {01..24}; do
    NanoPlot --fastq 05-fastq/barcode$i.fastq -o 02-qc/Nanoplot_barcode_raw$i -t 64
done

##################################################################################
# 4. Filter the fastq files for length and quality
mkdir -p 06-filtered
for i in {01..24}; do
    filtlong --min_length 1000 --min_mean_q 9 05-fastq/barcode$i.fastq > 06-filtered/barcode$i.fastq
done

# check fastq files for quality after filtering
for i in {01..24}; do
    echo "Stats for barcode$i.fastq" >> 02-qc/all_barcodes_afterfiltering.stats.txt
    seqkit stats 06-filtered/barcode$i.fastq >> 02-qc/all_barcodes_afterfiltering.stats.txt
    echo -e "\n" >> 02-qc/all_barcodes_afterfiltering.stats.txt # Add a newline for readability
done

# also check quality with NanoPlot
for i in {01..24}; do
    NanoPlot --fastq 06-filtered/barcode$i.fastq -o 02-qc/Nanoplot_barcode_filtered$i -t 64
done

##################################################################################
# 5. Assembly of the reads
# using flye to assemble the reads with nano-hq option
mkdir -p 07-assembly
for i in {01..24}; do
    flye --nano-hq 06-filtered/barcode$i.fastq --out-dir 07-assembly/barcode$i \
    -t 64 --iterations 3
done

##################################################################################
# some assemblies did not work, probably because of too many reads
# so i filtered more strictly and just use a part of the data

mkdir -p 06-filtered2
# barcode 03
# Leuconostoc: Genome size ~ 1.8Mb, total bases before filtering: ~ 1.5Gb, estimated coverage ~ 833x
# keep 1/8.3 of the data -> 12%
filtlong --min_length 1000 --min_mean_q 9 --keep_percent 12 05-fastq/barcode03.fastq > 06-filtered2/barcode03.fastq
# barcode 04
# Lactococcus: Genome size ~ 2.4Mb, total bases before filtering: ~ 2.3Gb, estimated coverage ~ 958x
# keep 1/9.6 of the data -> 10%
filtlong --min_length 1000 --min_mean_q 9 --keep_percent 10 05-fastq/barcode04.fastq > 06-filtered2/barcode04.fastq
# barcode 05
# Carnobacterium (divergens): Genome size ~ 3.2, total bases before filtering: ~ 2.6Gb, estimated coverage ~ 813x
# keep 1/8.1 of the data -> 12%
filtlong --min_length 1000 --min_mean_q 9 --keep_percent 12 05-fastq/barcode05.fastq > 06-filtered2/barcode05.fastq
# barcode 09
# Lactococcus: Genome size ~ 2.4Mb, total bases before filtering: ~ 2.0Gb, estimated coverage ~ 833x
# keep 1/8.3 of the data -> 12%
filtlong --min_length 1000 --min_mean_q 9 --keep_percent 12 05-fastq/barcode09.fastq > 06-filtered2/barcode09.fastq
# barcode 10
# Latilactobacillus: Genome size ~ 1.9Mb, total bases before filtering: ~ 795Mb, estimated coverage ~ 418x
# keep 1/4.2 of the data -> 24%
filtlong --min_length 1000 --min_mean_q 9 --keep_percent 24 05-fastq/barcode10.fastq > 06-filtered2/barcode10.fastq

# run the assembly of the additionally filtered data
for i in 03 04 05 09 10; do
    flye --nano-hq 06-filtered2/barcode$i.fastq --out-dir 07-assembly/barcode$i \
    -t 64 --iterations 3
done
conda deactivate

##################################################################################
# 6. Polishing the assembly
conda activate medaka

# using medaka for additional polishing
mkdir -p 08-polished
for i in 01 02 06 07 08 11 12 13 14 15 16; do
    medaka_consensus -i 06-filtered/barcode$i.fastq -d 07-assembly/barcode$i/assembly.fasta -o 08-polished/barcode$i -t 64 -m r1041_e82_400bps_sup_v5.0.0
done

for i in 03 04 05 09 10; do
    medaka_consensus -i 06-filtered2/barcode$i.fastq -d 07-assembly/barcode$i/assembly.fasta -o 08-polished/barcode$i -t 64 -m r1041_e82_400bps_sup_v5.0.0
done

mkdir 09-consensusfastas
for i in {01..16}; do
    if [ -f "08-polished/barcode$i/consensus.fasta" ]; then
        cp "08-polished/barcode$i/consensus.fasta" "09-consensusfastas/barcode$i.fasta"
        echo "Copied barcode$i/consensus.fasta to 09-consensusfastas/barcode$i.fasta"
    else
        echo "consensus.fasta not found in barcode$i/"
    fi
done

conda deactivate

##################################################################################
# 7. Assess the quality of the assemblies
conda activate Minion-2024
mkdir -p 10-checkm
checkm lineage_wf -t 64 09-consensusfastas/ 10-checkm/ -x fasta

checkm qa 10-checkm/lineage.ms 10-checkm/ > 10-checkm/quality_assessment.txt

conda deactivate

##################################################################################
# busco analysis of the assemblies
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/06-Minion
conda activate busco
mkdir -p 15-busco
for i in 03 04 09 12; do # leuconostoc, lactococcus, 
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l lactobacillaceae_odb12 -o 15-busco/barcode$i -c 100
done

for i in 10 14; do # latilactobacillus
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l lactobacillus_odb12 -o 15-busco/barcode$i -c 100
done

for i in 05 06; do # carnobacterium
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l carnobacterium_odb12 -o 15-busco/barcode$i -c 100
done

for i in 01; do # brochothrix
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l listeriaceae_odb12 -o 15-busco/barcode$i -c 100
done

for i in 02 11 07 08; do # Serratia, Kluyvera, Butauxiella
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l enterobacteriaceae_odb12 -o 15-busco/barcode$i -c 100
done

for i in 16; do # Pseudomonas
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l pseudomonas_odb12 -o 15-busco/barcode$i -c 100
done

for i in 15; do # Psychrobacter
    busco -i 09-consensusfastas/barcode${i}.fasta -m genome -l psychrobacter_odb12 -o 15-busco/barcode$i -c 100
done
conda deactivate
##################################################################################
# 8. Taxonomic classification of the assemblies using GTDBTk

# perparing a new environment
#conda create -n gtdbtk-2.4.0 -c conda-forge -c biodonda gtdbtk=2.4.0
#conda activate gtdbtk-2.4.0

# download the GTDB database to data
#download-db.sh /data/Unit_LMM/selberherr-group/roch/00_databases/gtdb_db_R220

# run GTDBTk
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/06-Minion
cd $workdir
mkdir -p 11-gtdbtk
gtdbtk classify_wf \
    --genome_dir 09-consensusfastas \
    --out_dir 11-gtdbtk \
    --cpus 64 \
    --mash_db /data/Unit_LMM/selberherr-group/roch/00_databases/mash_db/gtdbtk_database.msh \
    -x fasta


##################################################################################
# 9. using prokka to annotate the genomes

conda activate prokka-2024

mkdir 11-prokka

for genome in 09-consensusfastas/*.fasta; do
    sample=$(basename $genome .fasta)
    prokka \
    --outdir 11-prokka/${sample}_annotated \
    --prefix ${sample} \
    --kingdom Bacteria \
    --cpus 64 \
    --addgenes \
    --rfam \
    $genome
done

##################################################################################
# 10. using Antismash to identify the gene clusters
conda activate antismash-2024
mkdir -p 12-antismash

for i in {01..16}; do
  antismash --cb-general \
  --cb-knownclusters \
  --cb-subclusters \
  --asf \
  --pfam2go \
  --smcog-trees \
  --cpus 100 \
  --output-dir 12-antismash/barcode${i} \
  11-prokka/barcode${i}_annotated/barcode${i}.gbk
done

##################################################################################
# 11. using eggNOG-mapper to annotate the genomes
conda activate eggnog-2024
mkdir 13-eggnog

for i in {01..16}; do
    emapper.py -i 11-prokka/barcode${i}_annotated/barcode${i}.faa \
    --output barcode${i} \
    --output_dir 13-eggnog \
    --cpu 100 \
    --dmnd_db /data/Unit_LMM/selberherr-group/roch/00_databases/eggnog/eggnog_proteins.dmnd
done

# export the KEGG_ko to a txt
for i in {01..16}; do
    awk -F'\t' '!/^#/ && $12 != "-" && $12 != "" {
        gsub(/ko:/, "", $12);
        split($12, kos, ",");
        for (i in kos) {
            if (kos[i] ~ /^K/) {
                print $1 "\t" kos[i];
            }
        }
    }' 13-eggnog/barcode${i}.emapper.annotations | sort | uniq > 13-eggnog/kegg_ko_barcode${i}.txt
done

for i in {01..16}; do
    awk -F'\t' '!/^#/ && $11 != "-" && $11 != "" {
        split($11, ecs, ",");
        for (i in ecs) {
            print $1 "\t" ecs[i];
        }
    }' 13-eggnog/barcode${i}.emapper.annotations | sort | uniq > 13-eggnog/ec_barcode${i}.txt
done
##################################################################################
# 12. using abricate to identify the resistance genes
conda activate abricate
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/06-Minion
mkdir 14-resfinder
abricate 11-prokka/*/*.fna > 14-resfinder/resfinder_results_ncbi.tab \
    --db ncbi \
    --threads 64 

abricate 11-prokka/*/*.fna > 14-resfinder/resfinder_results_megares.tab \
    --db megares \
    --threads 64 


abricate 11-prokka/*/*.fna > 14-resfinder/resfinder_results_card.tab \
    --db card \
    --threads 64 

abricate 11-prokka/*/*.fna > 14-resfinder/resfinder_results_resfinder.tab \
    --db resfinder \
    --threads 64 
