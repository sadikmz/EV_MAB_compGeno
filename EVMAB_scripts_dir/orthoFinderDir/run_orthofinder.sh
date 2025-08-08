#!/bin/bash
#SBATCH --job-name=np_nwc1059
#SBATCH --partition=largemem
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
module load RAxML/8.2.12-gompi-2020a-hybrid-avx2

PROT_DIR="/home/data/salixomics/sadik/evmab_dir/ev_mab_prots_dmdb"
PATH_GENOMES="/home/data/salixomics/sadik/evmab_dir/genome_dir"
QPROT_DIR="/home/data/salixomics/sadik/evmab_dir/EV_MAB_novel_genes_db"
dmnd_blastdb=/home/data/salixomics/sadik//database_dir/NCBI_NR
threads=80

orthofinder -f prots_data -t $threads -a 8 -M msa -A muscle -T raxml

