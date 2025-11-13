#!/bin/bash
#SBATCH --job-name=mbev
#SBATCH --partition=compute
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

# module purge
#module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281
# module load picard/3.0.0-Java-17
module load Trim_Galore/0.6.10-GCCcore-11.3.0

#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREADS

WORKDIR=$(pwd)
# PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
#PICARD=/home/u1866313/miniconda3/envs/picard/share/picard-2.27.1-0
PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
query_genome="MB"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
mkdir "$PATH_WGS"
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps/
ref_genome=EV_bedadeti
cpus=28



#$APPSDIR/bwa-mem2/bwa-mem2 index ${PATH_GENOMES}/${ref_genome}.fna

for READS in ERR10695529 ERR10695603 ERR10695605 ERR10695608

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
	esearch -db sra -query PRJEB58004 | efetch -format runinfo | cut -d ',' -f 1 | grep -E $READS | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

	# Trim adapter sequences 

	trim_galore -cores 8 --quality 30  --output_dir $PATH_WGS --paired "$READS"_?.fastq.gz

	## extract coordinates of geneic and repeat regions 
	#gff2bed < $PATH_GENOMES/${ref_genome}.gff3 | awk '{if($8=="gene") print $1,$2,$3}' OFS='\t'  > $PATH_GENOMES/${ref_genome}.gene.bed
	#samtools faidx $ref
	#cut -f1,2 "$ref".fai  > "$ref".sizes 
	#bedtools complement -i $PATH_GENOMES/${ref_genome}.gene.bed -g "$ref".sizes > $PATH_GENOMES/${ref_genome}.gene.complement.bed 
	#cat $PATH_GENOMES/${ref_genome}.gene.bed $PATH_GENOMES/${ref_genome}.gene.complement.bed | bedtools sort > $PATH_GENOMES/${ref_genome}.gene_repeats.bed

	## map PE reads  
	# $APPSDIR/bwa-mem2/bwa-mem2 index $ref
	$APPSDIR/bwa-mem2-2.2.1_x64-linux/bwa-mem2  mem -t $cpus $ref $Fread $Rread | samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 | samtools sort - -@$cpus -o ${output_dir}/$out.allMapped.sorted.bam

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
	#java -Xmx100G -jar $PICARD/picard.jar MarkDuplicates \
	#	INPUT=${output_dir}/$out.allMapped.sorted.bam \
	#	O=${output_dir}/$out.allMapped.sorted.markdup.bam \
	#	M=${output_dir}/$out.allMapped.sorted.markdup.bam.metrics.txt \
	#	REMOVE_DUPLICATES=True

	rm ${output_dir}/$out.allMapped.sorted.bam

	## index bam files 
	samtools index ${output_dir}/$out.allMapped.sorted.markdup.bam -@ $cpus

	# qualimap

	module load Qualimap/2.2.1-foss-2021b-R-4.1.2

	qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.markdup.bam -outdir ${output_dir}/qualimap -outfile ${READS}_"$ref_genome".qualimap -sd -c -nt $cpus -outformat PDF:HTML -ip --java-mem-size=200G
	mv ${output_dir}/qualimap/genome_results.txt ${output_dir}/qualimap/${READS}_genome_results.txt 

	qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.markdup.bam -outdir ${output_dir}/qualimap_by_region -outfile ${READS}_"$ref_genome".qualimap_genes -sd -c -nt $cpus -gff $PATH_GENOMES/${ref_genome}.gene.bed -oc ${READS}_"$ref_genome".qualimap_genes_cov.txt -os ${READS}_"$ref_genome".qualimap_repeats -outformat PDF:HTML -ip --java-mem-size=200G
	mv ${output_dir}/qualimap_by_region/genome_results.txt ${output_dir}/qualimap_by_region/${READS}_genome_results.txt 


done 


./reads_mapping.01.sh
./reads_mapping.02.sh
./reads_mapping.03.sh
./reads_mapping.04.sh
./reads_mapping.05.sh
