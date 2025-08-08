#!/bin/bash
#SBATCH --job-name=mbev024
#SBATCH --partition=largemem
#SBATCH --time=4-00:00:00
#SBATCH --mem=2500G
#SBATCH --qos=generous-compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=79
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-user=sadik.abrare@rothamsted.ac.uk
#SBATCH --mail-type=END,FAIL

WORKDIR=$(pwd)
PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
query_genome="MB"
PATH_WGS="$WORKDIR"/"$query_genome"_trim_out
#mkdir "$PATH_WGS"
PATH_GENOMES=/home/data/salixomics/sadik/evmab_dir/genome_dir
APPSDIR=/home/data/salixomics/sadik/apps/
ref_genome=EV_mazia
cpus=79
ls -lht panMB_EV_mazia_out/*.bam |awk '{print $10}' > panMB_mazia.bamList.txt

#$APPSDIR/bwa-mem2/bwa-mem2 index ${PATH_GENOMES}/${ref_genome}.fna

for i in $(cat panMB_mazia.bamList.txt )

do
        READS=$(basename "$i" | sed 's/.allMapped.*//g' )	
	#query_genotype=$READS
	#mkdir $WORKDIR/panMB_${ref_genome}_out
	output_dir=$WORKDIR/panMB_${ref_genome}_out
	#Fread=${PATH_WGS}/${query_genotype}_1.fastq.gz
	#Rread=${PATH_WGS}/${query_genotype}_2.fastq.gz
	#ref=${PATH_GENOMES}/${ref_genome}.fna
	# cpus=48
	out=${ref_genome}_${READS}

	# download reads
	#esearch -db sra -query PRJEB58004 | efetch -format runinfo | cut -d ',' -f 1 | grep -E $READS | xargs -n 1 -P $cpus fastq-dump --split-files --gzip --skip-technical

	# Trim adapter sequences 

#	trim_galore -cores $cpus --quality 30  --output_dir $PATH_WGS --paired "$READS"_?.fastq.gz

	## extract coordinates of geneic and repeat regions 
	#gff2bed < $PATH_GENOMES/${ref_genome}.gff3 | awk '{if($8=="gene") print $1,$2,$3}' OFS='\t'  > $PATH_GENOMES/${ref_genome}.gene.bed
	#samtools faidx $ref
	#cut -f1,2 "$ref".fai  > "$ref".sizes 
	#bedtools complement -i $PATH_GENOMES/${ref_genome}.gene.bed -g "$ref".sizes > $PATH_GENOMES/${ref_genome}.gene.complement.bed 
	#cat $PATH_GENOMES/${ref_genome}.gene.bed $PATH_GENOMES/${ref_genome}.gene.complement.bed | bedtools sort > $PATH_GENOMES/${ref_genome}.gene_repeats.bed

	## map PE reads  
	# $APPSDIR/bwa-mem2/bwa-mem2 index $ref
	#$APPSDIR/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t $cpus $ref $Fread $Rread -P | samtools view - -Sb -@$cpus | samtools view -b -@$cpus -F 4 | samtools sort - -@$cpus -o ${output_dir}/$out.allMapped.sorted.bam

	# REMOVE_DUPLICATES 
        #samtools collate -o ${output_dir}/$out.allMapped.sorted.collate.bam ${output_dir}/$out.allMapped.sorted.bam -@$cpus 
	#rm ${output_dir}/$out.allMapped.sorted.bam

        #samtools fixmate -m ${output_dir}/$out.allMapped.sorted.collate.bam ${output_dir}/$out.allMapped.sorted.fixmate.bam -@$cpus
        #rm ${output_dir}/$out.allMapped.sorted.collate.bam

        #samtools sort -o ${output_dir}/$out.allMapped.sorted.02.bam ${output_dir}/$out.allMapped.sorted.fixmate.bam -@$cpus
        #rm ${output_dir}/$out.allMapped.sorted.fixmate.bam

        #samtools markdup ${output_dir}/$out.allMapped.sorted.02.bam ${output_dir}/$out.allMapped.sorted.markdup.bam -@$cpus
        #rm ${output_dir}/$out.allMapped.sorted.sorted.02.bam
	
	# java -Xmx100G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	#java -Xmx100G -jar $PICARD/picard.jar MarkDuplicates \
	#	INPUT=${output_dir}/$out.allMapped.sorted.bam \
	#	O=${output_dir}/$out.allMapped.sorted.markdup.bam \
	#	M=${output_dir}/$out.allMapped.sorted.markdup.bam.metrics.txt \
	#	REMOVE_DUPLICATES=True

	#rm ${output_dir}/$out.allMapped.sorted.bam

	## index bam files 
	#samtools index ${output_dir}/$out.allMapped.sorted.markdup.bam -@ $cpus

	# qualimap
	#module load Qualimap/2.2.1-foss-2021b-R-4.1.2

	#qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.markdup.bam -outdir ${output_dir}/qualimap -outfile ${READS}_"$ref_genome".qualimap -sd -c -nt $cpus -outformat PDF:HTML -ip --java-mem-size=200G
	#mv ${output_dir}/qualimap/genome_results.txt ${output_dir}/qualimap/${READS}_genome_results.txt 

	#qualimap bamqc -bam ${output_dir}/$out.allMapped.sorted.markdup.bam -outdir ${output_dir}/qualimap_by_region -outfile ${READS}_"$ref_genome".qualimap_genes -sd -c -nt $cpus -gff $PATH_GENOMES/${ref_genome}.gene.01.bed -oc ${READS}_"$ref_genome".qualimap_genes_cov.txt -os ${READS}_"$ref_genome".qualimap_repeats -outformat PDF:HTML -ip --java-mem-size=200G
	#mv ${output_dir}/qualimap_by_region/genome_results.txt ${output_dir}/qualimap_by_region/${READS}_genome_results.txt 

        ## Coverage all mapped reads in genic regions

        #bamCoverage -b ${output_dir}/$out.allMapped.sorted.markdup.bam --numberOfProcessors $cpus -o ${output_dir}/$out.allMapped.coverage.bw
	#echo $i >> test.txt
	#echo $READS >> test.txt
	bedtools bamtobed -i "$i"  | \
	bedtools coverage -a ${ref_genome}.genome_file.txt -iobuf 1000000M -b - > ${output_dir}/${READS}.allMapped.bed
done


