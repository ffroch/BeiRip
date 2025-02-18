# download genome files for bos taurus
cd /data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip/05-qiime
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCA_000003055.5_Bos_taurus_UMD_3.1.1/GCA_000003055.5_Bos_taurus_UMD_3.1.1_genomic.fna.gz
gzip -d GCA_000003055.5_Bos_taurus_UMD_3.1.1_genomic.fna.gz
conda activate blast
# make blast database
makeblastdb \
  -in GCA_000003055.5_Bos_taurus_UMD_3.1.1_genomic.fna \
  -out bos_taurus_DB \
  -dbtype nucl \
  -hash_index