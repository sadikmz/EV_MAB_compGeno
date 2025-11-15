#!/bin/bash
#mkdir fun_annot_out

# build local dbs 
makeblastdb in ~/data/uniprot/uniprot_sprot.fasta -dbtype prot
gunzip ~data/uniprot_trembl.fasta.gz 
makeblastdb -in ~data/uniprot_trembl.fasta -dbtype prot 
makeblastdb -in ~/data/cog_db/cog-20.fa -dbtype prot

# path to EggNOG diamond db
export EGGNOG_DATA_DIR=~/data/EggNOG
data_dir=~/data/
cpu=n
# Download eggNOG database
~/apps/eggnog/eggnog-mapper-2.1.8/download_eggnog_data.py -P -H -d 33090 -f --data_dir ~/data/EggNOG -y -q

# Run 
genotype="genotype_prefix"

# NCBI non-redundant db 

blastp \
-db ~/NCBI_NR/BLASTDB/nr \
-query "$genotype".25.TE_excluded.fasta \
-out "$genotype".maker2uni.cog.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads $cpu

# UNIPROT
blastp \
-db ~/data/uniprot/uniprot_sprot.fasta \
-query "$genotype".25.TE_excluded.fasta \
-out "$genotype".maker2uni.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads $cpu

# TeEMBL	
blastp \
-db ~data/uniprot_trembl.fasta \
-query "$genotype".25.TE_excluded.fasta \
-out "$genotype".maker2uni.trembl.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads $cpu

# cog

blastp \
-db ~/data/cog_db/cog-20.fa \
-query "$genotype".25.TE_excluded.fasta \
-out "$genotype".maker2uni.cog.blastp.txt \
-evalue 0.000001 \
-outfmt 6 \
-num_alignments 1 \
-seg yes \
-soft_masking true \
-lcase_masking \
-max_hsps 1 \
-num_threads $cpu


# Eggnog

emapper.py \
--cpu $cpu \
-m diamond \
-o eggNOG \
-i "$data_dir"/proteins/musa_ba.fasta \
--itype proteins \
--sensmode ultra-sensitive \
--output_dir eggnog_mapper_musa_ba.out  \
--dmnd_db "$data_dir"/EggNOG/eggnog_proteins.dmnd \
--evalue 1e-5 \
--override \
--report_orthologs 

