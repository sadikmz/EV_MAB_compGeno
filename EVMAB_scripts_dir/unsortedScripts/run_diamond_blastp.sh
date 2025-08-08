#!/bin/bash
#SBATCH --job-name=lastz
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
# ml BLAST+/2.14.1-gompi-2023a

PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0/
PROT_DIR="/home/data/salixomics/sadik/evmab_dir/EV_MAB_novel_genes_dmdb"
PATH_GENOMES="/home/data/salixomics/sadik/evmab_dir/genome_dir"
# ref_genome=EV_bedadeti
WORKDIR=$(pwd)
threads=32

/home/data/salixomics/sadik/apps/diamond makedb --in $PROT_DIR/EV_MAB.novel_genes.fa -d $PROT_DIR/EV_MAB.novel_genes.fa

/home/data/salixomics/sadik/apps/diamond blastp \
--db $PROT_DIR/EV_MAB.novel_genes.fa \
--query $PROT_DIR/EV_MAB.novel_genes.fa \
--ultra-sensitive \
--masking 0 \
--out EV_MAB.novel_genes.diamond.blastp.e05.usensitive.tsv \
--outfmt 6 qseqid qlen sallseqid slen qstart qend sstart send evalue bitscore length pident nident mismatch gapopen qcovhsp scovhsp \
--compress 1 \
--unal 1 \
--evalue 1e-5 \
--threads $threads 


# extract blast hits

zcat EV_MAB.novel_genes.diamond.blastp.e05.usensitive.tsv.gz | awk '{if ($5 > $6) print $1"\t"$6"\t"$5"\t"$12"\t"$11"\t"$10; else print $1"\t"$5"\t"$6"\t"$12"\t"$11"\t"$10}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' | bedtools sort | bedtools merge > EV_MAB.novel_genes.diamond.blastp.e05.usensitive.query.bed

zcat EV_MAB.novel_genes.diamond.blastp.e05.usensitive.tsv.gz | awk '{if ($7 > $8) print $1"\t"$8"\t"$7"\t"$12"\t"$11"\t"$10; else print $1"\t"$7"\t"$8"\t"$12"\t"$11"\t"$10}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' | bedtools sort | bedtools merge > EV_MAB.novel_genes.diamond.blastp.e05.usensitive.ref.bed

## whole proteom 

PROT_DIR="/home/data/salixomics/sadik/evmab_dir/ev_mab_prots_dmdb"
PATH_GENOMES="/home/data/salixomics/sadik/evmab_dir/genome_dir"
QPROT_DIR="/home/data/salixomics/sadik/evmab_dir/EV_MAB_novel_genes_db"

# ref_genome=EV_bedadeti
WORKDIR=$(pwd)
threads=32

/home/data/salixomics/sadik/apps/diamond makedb --in $PROT_DIR/EV_MAB.prot.fa -d $PROT_DIR/EV_MAB.prot.fa 

/home/data/salixomics/sadik/apps/diamond blastp \
--db $PROT_DIR/EV_MAB.prot.fa \
--query $QPROT_DIR/EV_MAB.novel_genes.fa \
--ultra-sensitive \
--masking 0 \
--out EV_MAB.novel_genes_vs_EV_MAB_prots.diamond.blastp.e05.usensitive.tsv \
--outfmt 6 qseqid qlen sallseqid slen qstart qend sstart send evalue bitscore length pident nident mismatch gapopen qcovhsp scovhsp \
--compress 1 \
--unal 1 \
--evalue 1e-5 \
--threads $threads 


# extract blast hits

zcat EV_MAB.novel_genes_vs_EV_MAB_prots.diamond.blastp.e05.usensitive.tsv.gz |  awk '{if ($5 > $6) print $1"\t"$6"\t"$5"\t"$12"\t"$11"\t"$10; else print $1"\t"$5"\t"$6"\t"$12"\t"$11"\t"$10}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' | bedtools sort | bedtools merge > EV_MAB.novel_genes_vs_EV_MAB_prots.diamond.blastp.e05.usensitive.query.bed


zcat EV_MAB.novel_genes_vs_EV_MAB_prots.diamond.blastp.e05.usensitive.tsv | awk '{if ($7 > $8) print $1"\t"$8"\t"$7"\t"$12"\t"$11"\t"$10; else print $1"\t"$7"\t"$8"\t"$12"\t"$11"\t"$10}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' | bedtools sort | bedtools merge > EV_MAB.novel_genes_vs_EV_MAB_prots.diamond.blastp.e05.usensitive.ref.bed


