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
ml BLAST+/2.14.1-gompi-2023a

PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0/
PROT_DIR="/home/data/salixomics/sadik/evmab_dir/EV_MAB_novel_genes_db"
PATH_GENOMES="/home/data/salixomics/sadik/evmab_dir/genome_dir"
# ref_genome=EV_bedadeti
WORKDIR=$(pwd)
threads=32

# create blastdb 

#makeblastdb -in $PROT_DIR/EV_MAB.novel_genes.fa -input_type fasta -dbtype prot
# makeblastdb -in $PATH_GENOMES/EV_bedadeti.fna -input_type fasta -dbtype nuc

# run blastp 

blastp \
-query $PROT_DIR/EV_MAB.novel_genes.fa \
-db $PROT_DIR/EV_MAB.novel_genes.fa \
-outfmt 6 \
-evalue 1e-5 \
-subject_besthit 1 \
-num_threads $threads \
-out EV_MAB_novel_prot_self_blastp.e05.out

# Filter blast hits > 25% pident
awk '{if ($7 > $8) print $1"\t"$8"\t"$7"\t"$3"\t"$4"\t"$12; else print $1"\t"$7"\t"$8"\t"$3"\t"$4"\t"$12}' EV_MAB_novel_prot_self_blastp.e05.out | awk '{if ($3> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }'> EV_MAB_novel_prot_self_blastp.query.bed 
awk '{if ($9 > $10) print $1"\t"$10"\t"$9"\t"$3"\t"$4"\t"$12; else print $1"\t"$9"\t"$10"\t"$3"\t"$4"\t"$12}' EV_MAB_novel_prot_self_blastp.e05.out | awk '{if ($3> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }'> EV_MAB_novel_prot_self_blastp.ref.bed 


