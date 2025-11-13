#!/bin/bash
#SBATCH --job-name=12qv
#SBATCH --partition=standard
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
#module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281
# module load picard/3.0.0-Java-17
# module load Trim_Galore/0.6.10-GCCcore-11.3.0
module load Trim_Galore/0.6.10-GCCcore-11.3.0

#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREADS

WORKDIR=$(pwd)
# PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
query_genome="EV_reads"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
#mkdir "$PATH_WGS"
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps/
ref_genome=EV_bedadeti
cpus=28

zcat EV_reads_trim_out/SRR4304971_1_val_1.fq.gz EV_reads_trim_out/SRR4304993_1_val_1.fq.gz > EV_reads_trim_out/Erpha20_1_val_1.fq.gz
zcat EV_reads_trim_out/SRR4304971_2_val_2.fq.gz EV_reads_trim_out/SRR4304993_2_val_2.fq.gz > EV_reads_trim_out/Erpha20_2_val_2.fq.gz

for READS in Erpha20

do
        query_genotype=$READS
        mkdir $WORKDIR/panEV_${ref_genome}_out
        output_dir=$WORKDIR/panEV_${ref_genome}_out
        Fread=${PATH_WGS}/${query_genotype}_1_val_1.fq.gz
        Rread=${PATH_WGS}/${query_genotype}_2_val_2.fq.gz
        ref=${PATH_GENOMES}/${ref_genome}.fna
        # cpus=48
        out=${ref_genome}_${query_genotype}

        # download reads
        #esearch -db sra -query PRJNA344540 | efetch -format runinfo | cut -d ',' -f 1 | grep 'SRR' | grep -E $READS | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

        # Trim adapter sequences

        #trim_galore -cores $cpus --quality 30  --output_dir $PATH_WGS --paired "$READS"_?.fastq.gz

	$APPSDIR/bwa-mem2/bwa-mem2 mem -t $cpus $ref $Fread $Rread | samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 | samtools sort - -@$cpus -o ${output_dir}/$out.allMapped.sorted.bam

        # REMOVE_DUPLICATES

        # java -Xmx100G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        java -Xmx100G -jar $PICARD/picard.jar MarkDuplicates \
        INPUT=${output_dir}/$out.allMapped.sorted.bam \
        O=${output_dir}/$out.allMapped.sorted.markdup.bam \
        M=${output_dir}/$out.allMapped.sorted.markdup.bam.metrics.txt \
        REMOVE_DUPLICATES=True

        rm ${output_dir}/$out.allMapped.sorted.bam

        ## index bam files
        samtools index ${output_dir}/$out.allMapped.sorted.markdup.bam -@ $cpus

        # coverage
        module load deepTools/3.5.2-foss-2022a
        bamCoverage -b ${output_dir}/$out.allMapped.sorted.markdup.bam --numberOfProcessors $cpus -o ${output_dir}/$out.allMapped.coverage.bw

done 

