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
# ml BLAST+/2.14.1-gompi-2023a

## whole proteom

PROT_DIR="/home/data/salixomics/sadik/evmab_dir/ev_mab_prots_dmdb"
PATH_GENOMES="/home/data/salixomics/sadik/evmab_dir/genome_dir"
QPROT_DIR="/home/data/salixomics/sadik/evmab_dir/EV_MAB_novel_genes_db"
dmnd_blastdb=/home/data/salixomics/sadik//database_dir/NCBI_NR
uniprot_dir="/home/data/salixomics/sadik/database_dir/uniprot_sprot/"
diamond="/home/data/salixomics/sadik/apps/diamond"
# ref_genome=EV_bedadeti
threads=80

# /home/data/salixomics/sadik/apps/diamond makedb --in $PROT_DIR/EV_MAB.prot.fa -d $PROT_DIR/EV_MAB.prot.fa
#/home/data/salixomics/sadik/apps/diamond prepdb --db "$dmnd_blastdb"/nr.gz --threads $threads 
"$diamond" makedb -in "$uniprot_dir"/athaliana_osativa/athalina_osativa.prot.fa --d "$uniprot_dir"/athaliana_osativa/athalina_osativa.prot.fa

"$diamond" blastp \
--db "$uniprot_dir"/athaliana_osativa/athalina_osativa.prot.fa \
--query $QPROT_DIR/EV_MAB.novel_genes.fa \
--ultra-sensitive \
--masking 0 \
--out EV_MAB.novel_genes.atos_db.blastp.e05.usensitive.tsv \
--outfmt 6 qseqid qlen sallseqid slen qstart qend sstart send evalue bitscore length pident nident mismatch gapopen qcovhsp scovhsp stitle \
--compress 1 \
--unal 1 \
--evalue 1e-5 \
--threads $threads

# extract blast hits
zcat EV_MAB.novel_genes.atos_db.blastp.e05.usensitive.tsv.gz | awk '{if ($5 > $6) print $1"\t"$6"\t"$5"\t"$12"\t"$11"\t"$10"\t"$11/$2"\t"$3"\t"$18; else print $1"\t"$5"\t"$6"\t"$12"\t"$11"\t"$10"\t"$11/$2"\t"$3"\t"$18}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  OFS='\t' > EV_MAB.novel_genes.atos_db.blastp.e05.usensitive.query.bed
zcat EV_MAB.novel_genes.atos_db.blastp.e05.usensitive.tsv.gz | awk '{if ($7 > $8) print $3"\t"$8"\t"$7"\t"$12"\t"$11"\t"$10"\t"$11/$4"\t"$1"\t"$18; else print $3"\t"$7"\t"$8"\t"$12"\t"$11"\t"$10"\t"$11/$4"\t"$1"\t"$18}' |  awk '{if ($4> 25) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  OFS='\t' > EV_MAB.novel_genes.atos_db.blastp.e05.usensitive.ref.bed

