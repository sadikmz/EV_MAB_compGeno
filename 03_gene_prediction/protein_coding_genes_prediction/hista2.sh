#!/bin/bash

genotype_prefix
export OMP_NUM_THREADS=n

RNAseq=~/data/rna_seq/DP8400009186BL_L01
cpus=$cpus
softmasked_genome=softmasked.fna # 
BASENAME=$(echo $softmasked_genome | sed "s/.fna//g")

## RNA-seq alignment using hisat2 (v2.1.0)

hisat2-build  $softmasked_genome ${softmasked_genome}.index -p $cpus --large-index

hisat2 \
-q -x ${softmasked_genome}.index \
-p $cpus \
-1 ${RNAseq}_573_1.fq.gz,${RNAseq}_574_1.fq.gz,${RNAseq}_575_1.fq.gz \
-2 ${RNAseq}_573_2.fq.gz,${RNAseq}_574_2.fq.gz,${RNAseq}_575_2.fq.gz \
--dta --dta-cufflinks \
--max-intronlen 160000 \
--no-mixed \
--very-sensitive \
--no-discordant \
--summary-file summary_file | samtools view -@$cpus -bS  | samtools sort -@$cpus -o ${BASENAME}.rnaseq.bam -T hisat2_tmp


stringtie -o ${BASENAME}_rnaseq.gff ${BASENAME}.rnaseq.bam  -p $cpus -A gene_abundance.out

