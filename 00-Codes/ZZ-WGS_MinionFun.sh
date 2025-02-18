
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/07-MinionFungi
conda activate Minion-2024
for i in {01..16};do 
    samtools fastq ./01-demux/*$i.bam > ./02-fastq/barcode$i.fq
done


cat barcode01.fq barcode02.fq barcode03.fq > sample0367.fq
cat barcode04.fq barcode05.fq > sample0991.fq
cat barcode06.fq barcode07.fq barcode08.fq > sample0368.fq
cat barcode09.fq barcode10.fq > sample1109.fq
cat barcode11.fq barcode12.fq > sample0317.fq
cat barcode13.fq barcode14.fq barcode15.fq barcode15.fq > sample0049.fq

mkdir 03-qc
# check fastq files for quality
for i in 0367 0991 0368 1109 0317 0049; do
    echo "Stats for sample$i.fq" >> 03-qc/all_samples.stats.txt
    seqkit stats 02-fastq/sample$i.fq >> 03-qc/all_samples.stats.txt
    echo -e "\n" >> 03-qc/all_samples.stats.txt # Add a newline for readability
done

# filter the fastq files for length and quality
mkdir -p 04-filtered
for i in 0367 0991 0368 1109 0317; do
    filtlong --min_length 1000 --min_mean_q 9 02-fastq/sample$i.fq > 04-filtered/sample$i.fastq
done

# expected genome size for 0367 (Kurtzmaniella) 15Mb, total bases after filtering: 173 Mb -> 11x coverage
# expected genome size for 0991 (Debaryomyces) 12Mb, total bases after filtering: 324 Mb -> 27x coverage
# expected genome size for 0368 (Debaryomyces) 12Mb, total bases after filtering: 284 Mb -> 24x coverage
# expected genome size for 1109 (Barnettozyma) 11Mb, total bases after filtering: 166 Mb -> 15x coverage
# expected genome size for 0317 (Yarrowia) 20Mb, total bases after filtering: 156 Mb -> 8x coverage

# Lactococcus: Genome size ~ 2.4Mb, total bases before filtering: ~ 2.6Gb, estimated coverage ~ 1080x
# keep 1/10.8 of the data -> 9%
for i in 0049; do
    seqkit rmdup -s -i -o 02-fastq/sample$i_cleaned.fq 02-fastq/sample$i.fq
    filtlong --min_length 1000 --min_mean_q 9 --keep_percent 9 02-fastq/sample$i_cleaned.fq > 04-filtered/sample$i.fastq
done


# check fastq files for quality after filtering
for i in 0367 0991 0368 1109 0317 0049; do
    echo "Stats for sample$i.fastq" >> 03-qc/all_samples_afterfiltering.stats.txt
    seqkit stats 04-filtered/sample$i.fastq >> 03-qc/all_samples_afterfiltering.stats.txt
    echo -e "\n" >> 03-qc/all_samples_afterfiltering.stats.txt # Add a newline for readability
done

# proceed with the yeasts only

# Step 1: Flye Assembly
mkdir -p 05-assembly
for i in 0367 0991 0368 1109 0317; do
    flye --nano-hq 04-filtered/sample$i.fastq --out-dir 05-assembly/sample$i --min-overlap 1000 \
    -t 100 --iterations 3
done

conda deactivate

conda activate medaka

# Step 2: Racon Polishing
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/07-MinionFungi
for i in 0367 0991 0368 1109 0317; do
    mkdir -p ${workdir}/06-polished2/sample$i
        contigs="${workdir}/05-assembly/sample${i}/assembly.fasta"
        threads=100
        reads="${workdir}/04-filtered/sample${i}.fastq"

        minimap2 -x ava-ont -t $threads $contigs $reads > ${workdir}/05-assembly/sample${i}/overlaps_1.paf
        racon -t $threads $reads ${workdir}/05-assembly/sample${i}/overlaps_1.paf $contigs > ${workdir}/05-assembly/sample${i}/racon_1.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/05-assembly/sample${i}/racon_1.fasta $reads > ${workdir}/05-assembly/sample${i}/overlaps_2.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/05-assembly/sample${i}/overlaps_2.paf ${workdir}/05-assembly/sample${i}/racon_1.fasta > ${workdir}/05-assembly/sample${i}/racon_2.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/05-assembly/sample${i}/racon_2.fasta $reads > ${workdir}/05-assembly/sample${i}/overlaps_3.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/05-assembly/sample${i}/overlaps_3.paf ${workdir}/05-assembly/sample${i}/racon_2.fasta > ${workdir}/05-assembly/sample${i}/racon_3.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/05-assembly/sample${i}/racon_3.fasta $reads > ${workdir}/05-assembly/sample${i}/overlaps_4.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/05-assembly/sample${i}/overlaps_4.paf ${workdir}/05-assembly/sample${i}/racon_3.fasta > ${workdir}/05-assembly/sample${i}/racon_4.fasta

        medaka_consensus -i $reads -d ${workdir}/05-assembly/sample${i}/racon_4.fasta -o ${workdir}/06-polished2/sample${i} -t $threads -m r1041_e82_400bps_sup_v5.0.0
        rm -rf ${workdir}/05-assembly/sample${i}/*.paf ${workdir}/05-assembly/sample${i}/racon_*.fasta*
done


conda deactivate

conda activate busco
mkdir -p 07-busco
for i in 0367 0991 0368 1109 0317; do
    busco -i 06-polished2/sample$i/consensus.fasta -m genome -l ascomycota_odb10 -o 07-busco/sample$i -c 100
done

# check if the missing buscos from the both debaryomyces are the same
comm -12 >(sort 07-busco/sample0368/run_ascomycota_odb10/missing_busco_list.tsv) <(sort 07-busco/sample0991/run_ascomycota_odb10/missing_busco_list.tsv)
# 20 of 21 missing buscos are the same
# check if the missing buscos from Kurtzmaniella and Barnettozyma are the same
comm -12 >(sort 07-busco/sample1109/run_ascomycota_odb10/missing_busco_list.tsv) <(sort 07-busco/sample0367/run_ascomycota_odb10/missing_busco_list.tsv)
# 67 of 68/108 missing buscos are the same
# and now a cross check between 0368 and 0367
comm -12 <(sort 07-busco/sample0368/run_ascomycota_odb10/missing_busco_list.tsv) <(sort 07-busco/sample0367/run_ascomycota_odb10/missing_busco_list.tsv)
# 11 of the missing buscos are the same

mkdir -p 07-busco2
for i in 0367 0991 0368 1109 0317; do
    busco -i 06-polished2/sample$i/consensus.fasta -m genome -l saccharomycetes_odb10 -o 07-busco2/sample$i -c 100
done

conda deactivate

## before proceeding use kraken2 for taxonomic classification
conda activate kraken2
mkdir -p 08-kraken2
for i in 0367 0991 0368 1109; do
    kraken2 --db /data/Unit_LMM/databases/kraken2/k2_pluspf_20240409 --threads 100 --report 08-kraken2/sample${i}.report --output 08-kraken2/sample${i}.kraken2.out 06-polished2/sample$i/consensus.fasta
done
conda deactivate
#######################################################################
# run R code to combine the kraken2 reports with flye assembly stats
# run R code ZZ-WGSAssemblyCleaning.R
#######################################################################
# the fasta assemblies are now filtered for unwanted species contigs, repeat busco to check for completeness
conda activate busco
mkdir -p 07-busco3
for i in 0367 0991 0368 1109; do
    busco -i 06-polished2/sample$i/consensus_filtered.fasta -m genome -l saccharomycetes_odb10 -o 07-busco3/sample$i -c 100
done
conda deactivate

mkdir -p 09-blast
for i in 0367 0991 0368 1109; do
    blastn -query 06-polished2/sample$i/consensus_filtered.fasta -db /data/Unit_LMM/selberherr-group/roch/00_databases/ncbi_fungi/fungi_db/fungi_db -out 09-blast/sample$i.blastn -num_threads 100 -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
done

cd ./09-blast
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz

for i in 0367 0991 0368 1109; do
    cut -f2 sample$i.blastn | sort | uniq > sample$i.unique_sseqids.txt
    grep -Fwf sample$i.unique_sseqids.txt nucl_gb.accession2taxid > sample$i.matched_taxids.txt
done

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz

awk -F '|' '$4 ~ /scientific name/ {
    gsub(/^[ \t]+|[ \t]+$/, "", $1);  # Trim leading/trailing spaces from taxid
    gsub(/^[ \t]+|[ \t]+$/, "", $2);  # Trim leading/trailing spaces from name
    print $1 "\t" $2
}' OFS='\t' names.dmp > taxid_to_name.txt

# Run R code
# "bringblastandtaxontogether.R"

# run annotation with companion: https://companion.gla.ac.uk/jobs/new
# ID-0367:
    # Job name: ID-0367
    # Species prefix: ID0367
    # Species name: [Candida] zeylanoides
    # Target sequence: 06-polished2/sample0367/consensus_filtered.fasta
    # Transcript evidence: No, do not use transcript evidence
    # Reference organism: Debaryomyces hansenii CBS767
    # Pseudochromosome contiguation: Yes, contiguate pseudochromosomes
        # minimum required match length for contig placement: 500bp
        # minimum required matdch similarity for contig placement: 85%
    # Advanced settings: default

# ID-0368:
    # Job name: ID-0368
    # Species prefix: ID0368
    # Species name: Debaryomyces hanesenii
    # Target sequence: 06-polished2/sample0368/consensus_filtered.fasta
    # Transcript evidence: No, do not use transcript evidence
    # Reference organism: Debaryomyces hansenii CBS767
    # Pseudochromosome contiguation: Yes, contiguate pseudochromosomes
        # minimum required match length for contig placement: 500bp
        # minimum required matdch similarity for contig placement: 85%
    # Advanced settings: default

# ID-0991:
    # Job name: ID-0991
    # Species prefix: ID0991
    # Species name: Debaryomyces hanesenii
    # Target sequence: 06-polished2/sample0991/consensus_filtered.fasta
    # Transcript evidence: No, do not use transcript evidence
    # Reference organism: Debaryomyces hansenii CBS767
    # Pseudochromosome contiguation: Yes, contiguate pseudochromosomes
        # minimum required match length for contig placement: 500bp
        # minimum required matdch similarity for contig placement: 85%
    # Advanced settings: default

# ID-1109:
    # Job name: ID-1109
    # Species prefix: ID1109
    # Species name: [Candida] norvegica
    # Target sequence: 06-polished2/sample1109/consensus_filtered.fasta
    # Transcript evidence: No, do not use transcript evidence
    # Reference organism: Dyberlindnera fabianii
    # Pseudochromosome contiguation: Yes, contiguate pseudochromosomes
        # minimum required match length for contig placement: 500bp
        # minimum required matdch similarity for contig placement: 85%
    # Advanced settings: default

# Use proteins.fasta from companion to build model on modelseed.org
# "UPLOAD Microbes FASTA"
# Select Template Model: "Automatically select"
# Genome Type: "Microbial feature or Protein sequences"
# "Selected media: complete"

# I did the same for the bacterial faa files from prokka output


# sbml files were used in different combination in PyCoMo
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
for i in 0367 0368 0991 1109; do
    for j in 0049 0052 0136 0184 0467 0486 1047 1609 1626; do
    mkdir -p 11-sbml_models/${i}_${j} 
        pycomo --input 11-sbml_models/*${i}.sbml 11-sbml_models/*${j}.sbml \
                --name ${i}_${j} \
                --num-cores 100 \
                --output-dir 11-sbml_models/${i}_${j} \
                --fva-flux 0.9 \
                --fva-interaction \
                --composition-agnostic
    done
done


# sbml files were used in differenc combination with SMETANA
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip
mkdir -p 12-smetana
smetana 11-sbml_models/*0367.xml 11-sbml_models/*0052.xml \
    --output 12-smetana/0367_0052 \
    --flavor cobra \
    --detailed \
    --molweight \
    --solver cplex