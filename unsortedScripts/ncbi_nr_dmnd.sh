#!/bin/bash
#SBATCH --job-name=unip
#SBATCH --partition=largemem
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=generous-memory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

module load GCC/11.3.0 OpenMPI/4.1.4

# set variables 
threads=80
# uniprot_sprot=/home/data/salixomics/sadik/database_dir/uniprot_sprot/athaliana_osativa/athalina_osativa.prot.fa
#uniprot_sprot=/home/data/salixomics/sadik/database_dir/uniprot_sprot/uniprot_sprot.fasta.gz
ncbi_dmnd_nr=/home/data/salixomics/sadik/database_dir/NCBI_NR/nr.gz

# index database 
#diamond blastp diamond makedb --in $uniprot_sprot --db  $uniprot_sprot
diamond diamond prepdp -d $ncbi_dmnd_nr

# ~/apps/diamond blastp diamond makedb --in ~/data/uniprot/uniprot_sprot.fasta -d ~/data/uniprot/uniprot_sprot.fasta
# ~/apps/diamond diamond prepdp -d ~/data/NCBI_NR/nr.gz
predicted_proteins=/home/data/salixomics/sadik/evmab_dir/EV_MAB.nogo_noreads.prot.fasta
genotype=EV_MAB.nogo_noreads
#out=${genotype}.uniprot_sprot.diamond.blastp.e05.out 
out=${genotype}.ncbi_nr.diamond.blastp.e05.out 

# run blast
diamond blastp \
--query $predicted_proteins \
--db $ncbi_dmnd_nr \
--ultra-sensitive \
--masking none \
--strand both \
--out $out \
--outfmt 6 qseqid sseqid qlen slen qstart qend sstart send qseq sseq stitle pident approx_pident nident qstrand length evalue bitscore \
--header simple \
--compress 1 \
--unal 1 \
--evalue 1e-5 \
--threads $threads \
--no-parse-seqids
