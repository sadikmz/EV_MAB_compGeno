#!/bin/bash
#SBATCH --job-name=mzbdbd
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

#module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281
# module load picard/3.0.0-Java-17
# module load Trim_Galore/0.6.10-GCCcore-11.3.0
#module load Trim_Galore/0.6.10-GCCcore-11.3.0
#module load SAMtools/1.18-GCC-12.3.0
#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREADS

#module load Miniconda3/23.9.0-0
#module load Trim_Galore/0.6.10-GCCcore-11.3.0

#Activate base environmetn
#source ~/.bashrc

# Load bash_profile to avail to tools appeneded in the $PATH.
# Check the compatibilty of tools in bash_profile $PATH with those to be used in conda env before loading
#source ~/.bash_profile

# activate the specific env to run
#source activate /home/data/salixomics/sadik/miniconda3/envs/genome_qc

WORKDIR=$(pwd)
# PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
PICARD="/home/data/salixomics/sadik/miniconda3/envs/genome_qc/share/picard-3.2.0-0"
query_genome="EV_reads"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
mkdir "$PATH_WGS"
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps
#ref_genome=EV_bedadeti
cpus=79
genotype=bedadeti

#zcat EV_reads_trim_out/SRR151526?_2_val_2.fq.gz | bgzip -@$cpus > "$PATH_WGS"/"$genotype"_2_val_2.fq.gz  
#zcat EV_reads_trim_out/SRR151526?_1_val_1.fq.gz | bgzip -@$cpus > "$PATH_WGS"/"$genotype"_1_val_1.fq.gz

READS="Erpha20"
#for READS in SRR1515269	SRR1515268

for ref_genome in EV_mazia EV_bedadeti
do
        query_genotype=$READS
        #mkdir $WORKDIR/panEV_${ref_genome}_out
        output_dir=$WORKDIR/panEV_${ref_genome}_out
        Fread=${PATH_WGS}/${query_genotype}_1_val_1.fq.gz
        Rread=${PATH_WGS}/${query_genotype}_2_val_2.fq.gz
        ref=${PATH_GENOMES}/${ref_genome}.fna
        # cpus=48
        out=${ref_genome}_${query_genotype}

        # download reads
       # esearch -db sra -query PRJNA252658  | efetch -format runinfo | cut -d ',' -f 1 | grep -E $READS | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

        # Trim adapter sequences

        #trim_galore -cores $cpus --quality 30  --output_dir $PATH_WGS --paired "$READS"_?.fastq.gz
        ## map PE reads  
        # bwa-mem2 index $ref
        #$APPSDIR/bwa-mem2-2.2.1/bwa-mem2  mem -t $cpus $ref $Fread $Rread | samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 | samtools sort - -@$cpus -o ${output_dir}/$out.allMapped.sorted.bam

        # REMOVE_DUPLICATES 
        # REMOVE_DUPLICATES
        samtools collate -o ${output_dir}/$out.allMapped.sorted.collate.bam ${output_dir}/$out.allMapped.sorted.bam -@$cpus
        rm ${output_dir}/$out.allMapped.sorted.bam

        samtools fixmate -m ${output_dir}/$out.allMapped.sorted.collate.bam ${output_dir}/$out.allMapped.sorted.fixmate.bam -@$cpus
        rm ${output_dir}/$out.allMapped.sorted.collate.bam

        samtools sort -o ${output_dir}/$out.allMapped.sorted.02.bam ${output_dir}/$out.allMapped.sorted.fixmate.bam -@$cpus
        rm ${output_dir}/$out.allMapped.sorted.fixmate.bam

        samtools markdup ${output_dir}/$out.allMapped.sorted.02.bam ${output_dir}/$out.allMapped.sorted.markdup.bam -@$cpus
        rm ${output_dir}/$out.allMapped.sorted.sorted.02.bam

        # java -Xmx100G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
       # java -Xmx100G -jar $PICARD/picard.jar MarkDuplicates \
       # INPUT=${output_dir}/$out.allMapped.sorted.bam \
       # O=${output_dir}/$out.allMapped.sorted.markdup.bam \
       # M=${output_dir}/$out.allMapped.sorted.markdup.bam.metrics.txt \
       # REMOVE_DUPLICATES=True

#        rm ${output_dir}/$out.allMapped.sorted.bam

        ## index bam files 
        samtools index ${output_dir}/$out.allMapped.sorted.markdup.bam -@ $cpus

        # coverage 
        module load deepTools/3.5.2-foss-2022a
        bamCoverage -b ${output_dir}/$out.allMapped.sorted.markdup.bam --numberOfProcessors $cpus -o ${output_dir}/$out.allMapped.coverage.bw

done 

