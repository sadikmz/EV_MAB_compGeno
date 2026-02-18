# Genome QC and Profiling

This module covers genome size estimation, heterozygosity estimation, ploidy assessment, and assembly quality control for the enset (*Ensete ventricosum*) landraces Mazia and Bedadeti.

## Overview

Steps performed:

1. Prepare genomic reads (stLFR/NGS)
2. Estimate genome size and heterozygosity (jellyfish, GenomeScope2, FindGSE)
3. Assess ploidy level (Smudgeplot via FastK and PloidyPlot)
4. Assess assembly contiguity and gene-space completeness (gfastats, BUSCO)

## Requirements

- `jellyfish` >= 2.3.0
- `GenomeScope2`
- `FindGSE` (R package)
- `FastK`, `PloidyPlot`, `smudgeplot`
- `Trim Galore`
- `gfastats`
- `BUSCO` >= 5.2.1
- `NCBI Entrez Direct` (`esearch`, `efetch`, `fastq-dump`)
- `stLFRdenovo` (v1.0.5) for Mazia stLFR reads

Install conda-managed tools via the root [`environment.yml`](../environment.yml):

```bash
conda env create -f ../environment.yml
conda activate ev_mab_compgeno
```

---

### 1. Preparing genomic reads (stLFR/NGS)

```bash
#!/bin/bash

prog=~/apps/jellyfish/
wgs_dir=/home/data/wgs

# Download SRA files for landrace Bedadeti (SRR1515268, SRR1515269) and convert to FASTQ format
esearch -db sra -query PRJNA432894 \
  | efetch -format runinfo \
  | cut -d ',' -f 1 \
  | grep SRR \
  | grep -E 'SRR1515268|SRR1515269' \
  | xargs -n 1 -P 12 fastq-dump --split-files --gzip --skip-technical

zcat $wgs_dir/trim_out/SRR1515269_1.fastq.gz $wgs_dir/trim_out/SRR1515268_1.fastq.gz \
  > $wgs_dir/bedadeti_1.fastq.gz
zcat $wgs_dir/trim_out/SRR1515269_2.fastq.gz $wgs_dir/trim_out/SRR1515268_2.fastq.gz \
  > $wgs_dir/bedadeti_2.fastq.gz

rm $wgs_dir/trim_out/SRR151526?_?.fastq.gz

# For landrace Mazia, remove stLFR adapter sequences using
# stLFR2Supernova (v2.1.1): https://github.com/BGI-Qingdao/stlfr2supernova_pipeline
# and merge them with Mazia NGS reads

# Create output directory
mkdir stLFR_assemble
cd stLFR_assemble

# Set variables
data_dir=/home/data
out_path=output
max_memory=252G
threads=64
apps_dir=/home/apps/stLFRdenovo-v1.0.5

# Run stLFRdenovo
$apps_dir/stLFRdenovo_pipeline.sh \
  -f $data_dir/stLFR_raw_data/CL200149738_L01_read_1.fq.gz \
     $data_dir/stLFR_raw_data/CL200149738_L01_read_2.fq.gz \
  -o $out_path \
  -m $max_memory \
  -t $threads

# Combine adapter-trimmed stLFR reads with 100-bp raw paired-end reads
zcat stLFR_assemble/output/read.1.fq.gz.clean.gz \
     $data_dir/WGS_data/DP8400009105TR_L01_520_1.fq.gz \
  > $wgs_dir/mazia.1.fq.gz
zcat stLFR_assemble/output/read.2.fq.gz.clean.gz \
     $data_dir/WGS_data/DP8400009105TR_L01_520_2.fq.gz \
  > $wgs_dir/mazia.2.fq.gz

ls -lht $wgs_dir/*.?.fq.gz \
  | awk '{print $9}' \
  | sed 's/...fq.gz//g' \
  | uniq > accession_list

# accession_list should contain:
# mazia
# bedadeti
```

---

### 2. Genome size and heterozygosity estimation

```bash
#!/bin/bash
data_dir=/home/data

# Run jellyfish and GenomeScope2 across k-mer sizes of 17, 21, 27, and 31.
# Estimated genome size ~580 Mb, estimated coverage ~43x.
for K in 17 21 27 31; do
  for L in $(cat accession_list); do

    # Trim adapters and low-quality reads
    trim_galore \
      --cores 40 \
      --quality 30 \
      --fastqc \
      -o "${data_dir}_trim_out" \
      --paired \
      "${data_dir}/${L}.1.fq.gz" \
      "${data_dir}/${L}.2.fq.gz"

    # Generate k-mer table
    ${prog}/jellyfish-linux count \
      -C -m ${K} -s 580M --bf-size 25G -t 16 \
      <(zcat "${L}_1_val_1.fq.gz") \
      <(zcat "${L}_2_val_2.fq.gz") \
      -o "${L}.${K}_mer.reads.jf"

    # Generate k-mer histogram
    ${prog}/jellyfish-linux histo -t 16 \
      "${L}.${K}_mer.reads.jf" > "${L}.${K}_mer.histo"
    rm -rf "${L}.${K}_mer.reads.jf"

    # Run GenomeScope2
    genomescope2 \
      -i "${L}.${K}_mer.histo" \
      -o "genomescope_${L}_k${K}" \
      -p 2 \
      -k ${K} \
      --testing \
      --fitted_hist
  done
done

# Genome size estimation using FindGSE.
# K-mer histogram files generated above are used as input to the FindGSE R package.

Rscript << 'EOF'
library(FindGSE)
accessions <- c("mazia", "bedadeti")
kmer_sizes <- c(17, 21, 27, 31)
for (acc in accessions) {
  for (k in kmer_sizes) {
    histo_file <- paste0(acc, ".", k, "_mer.histo")
    output_prefix <- paste0("FindGSE_", acc, "_k", k)
    result <- FindGSE(histo_file, k = k, plot = TRUE)
    ggsave(filename = paste0(output_prefix, "_FindGSE_plot.png"), plot = last_plot())
    write.csv(result, file = paste0(output_prefix, "_FindGSE_results.csv"))
  }
}
EOF
```

---

### 3. Ploidy level estimation using Smudgeplot

```bash
#!/bin/bash

data_dir=/home/data/
tmp="/dev/shm"
threads=16
K=31

# Set to the filename prefix of trimmed reads for each genome
genotype=prefix_of_trimmed_reads_name_for_each_genome

# Create output directory
mkdir "${genotype}_output" "tmp_k${K}"

# Run FastK to create a k-mer database
FastK \
  -v -t4 -k${K} -M16 -T${threads} \
  "${data_dir}/trim_out/"*.gz \
  -N "${genotype}_output/FastK_Table" \
  -P "tmp_k${K}"

# Run PloidyPlot to find all k-mer pairs in the dataset.
# Generates <genotype>_output/kmerpairs_text.smu:
# a flat file with three columns: covB, covA, freq
PloidyPlot \
  -e12 -k -v -T${threads} \
  -o "${genotype}_output/kmerpairs" \
  "${genotype}_output/FastK_Table" \
  -P "${tmp}"

# Use the .smu file to infer ploidy and create Smudgeplot
smudgeplot.py plot \
  -n 15 \
  -t "${genotype}" \
  -o "${genotype}_output/${genotype}_smplot" \
  "${genotype}_output/kmerpairs_text.smu"

# Verify that 5 output files are generated (2 PDFs, 1 TSV summary, 2 TXT logs)
ls "${genotype}_output/${genotype}_smplot"*

# Remove the temporary directory created by FastK
rm -rf "tmp_k${K}"
```

---

### 4. Assembly QC

```bash
#!/bin/bash

# Input files:
#   mazia.fna
#   bedadeti.fna

for i in *.fna; do

  BASENAME=$(basename "${i}" .fna)

  # Assess assembly contiguity with gfastats
  ~/apps/gfastats/build/bin/gfastats \
    -b s \
    -f "${i}" \
    -s s \
    -t \
    --seq-report \
    --stats \
    --sort descending \
    > "${BASENAME}.gfastats.summary.txt"

  # Assess BUSCO gene-space completeness (embryophyta lineage)
  busco \
    -i "${i}" \
    -m genome \
    -l embryophyta_odb10 \
    -c 16 \
    -o "busco.${BASENAME}.out" \
    --long \
    -f

done
```
