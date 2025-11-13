#!/bin/bash
#SBATCH --job-name=mzmbpb
#SBATCH --partition=largemem
#SBATCH --time=7-00:00:00
#SBATCH --mem=2500G
#SBATCH --qos=generous-compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=79
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
#module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281
#module load picard/3.0.0-Java-17
#module load SAMtools/1.18-GCC-12.3.0
#module load Qualimap/2.2.1-foss-2021b-R-4.1.2
#module load Trim_Galore/0.6.10-GCCcore-11.3.0
module load minimap2/2.26-GCCcore-12.3.0


#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREADS

WORKDIR=$(pwd)
# PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
query_genome="MB"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps/
ref_genome=EV_mazia
cpus=79
threads=79
#zcat SRR24350321.fastq.gz SRR24235289.fastq.gz | bgzip -@$cpus > MA.pbhifi.fastq.gz

for READS in MB_pacbio
do 
query_genotype=$READS
#mkdir $WORKDIR/panEV_${ref_genome}_out
output_dir=$WORKDIR/panMB_${ref_genome}_out
Fread=${PATH_WGS}/${query_genotype}.1_val_1.fq.gz
Rread=${PATH_WGS}/${query_genotype}.2_val_2.fq.gz
hifi_reads=${WORKDIR}/bananaB.subreads.fastq.gz 
ref=${PATH_GENOMES}/${ref_genome}.fna
# cpus=48
out=${ref_genome}_${query_genotype}

#esearch -db sra -query PRJEB58004 | efetch -format runinfo | cut -d ',' -f 1 | grep ERR | grep -E $READS | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical


#cutadapt \
#-b "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT;min_overlap=35" \
#-b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=35" \
#--discard-trimmed \
#--revcomp \
#-e 0.1 \
#--report=minimal \
#--discard  \
#--output ${PATH_WGS}/"$query_genotype".hifi.adapt_trimmed.fastq.gz  \
#--cores ${threads} \
#"$query_genotype".fastq.gz


# minimap2 alignment 
minimap2 -ax map-pb ${PATH_GENOMES}/EV_mazia.fna $hifi_reads -t "$threads" -k 19 | samtools view - -Sb -@"$threads" | samtools sort -@"$threads" -o ${output_dir}/$out.allMapped.sorted.bam
#

samtools index ${output_dir}/$out.allMapped.sorted.bam

qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.bam -outdir ${output_dir}/qualimap -outfile ${READS}_"$ref_genome".qualimap -sd -c -nt $cpus -outformat PDF:HTML -ip --java-mem-size=200G
mv ${output_dir}/qualimap/genome_results.txt ${output_dir}/qualimap/${READS}_genome_results.txt 

qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.bam -outdir ${output_dir}/qualimap_by_region -outfile ${READS}_"$ref_genome".qualimap_genes -sd -c -nt $cpus -gff $PATH_GENOMES/${ref_genome}.gene.bed -oc ${READS}_"$ref_genome".qualimap_genes_cov.txt -os ${READS}_"$ref_genome".qualimap_repeats -outformat PDF:HTML -ip --java-mem-size=200G
mv ${output_dir}/qualimap_by_region/genome_results.txt ${output_dir}/qualimap_by_region/${READS}_genome_results.txt 

 ## Coverage all mapped reads in genic regions
bamCoverage -b ${output_dir}/$out.allMapped.sorted.bam --numberOfProcessors $cpus -o ${output_dir}/$out.allMapped.coverage.bw
done 
