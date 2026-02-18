# Comparative Genomics: Presence/Absence Variation and Species-Specific Gene Identification

Identification of genes present in *Ensete ventricosum* but absent from *Musa* species (and vice versa) using whole-genome read mapping and CDS alignment.

## Overview

The pipeline proceeds in three steps:

1. **Data preparation** - Extract genic and intergenic coordinates from genome annotations; index reference genomes for read mapping.
2. **Read mapping and PAV analysis** - Map whole-genome resequencing reads from multiple enset genotypes to reference genomes, assess coverage, and determine gene presence/absence.
3. **Identify species-specific genes** - Combine PAV evidence and CDS alignment results to define genes unique to enset or *Musa* species.

## Requirements

- `Trim Galore`
- `BWA-MEM2` >= 2.2.1
- `SAMtools` >= 1.22
- `Picard` >= 2.27.1
- `Qualimap` >= 2.2.2
- `BEDtools` >= 2.30.0
- `deepTools`
- `Python` >= 3.9 with `pybedtools`, `pysam`, `pandas`
- NCBI SRA tools (`fastq-dump`)

Install conda-managed tools via the root [`environment.yml`](../environment.yml):

```bash
conda env create -f ../environment.yml
conda activate ev_mab_compgeno
```

## Steps

### 1. Data preparation

[`00_prep_data.sh`](00_prep_data.sh) - Extract genic and intergenic coordinates from genome annotations, and index reference genomes for read mapping.

### 2. Read mapping and PAV analysis

[`01_readsCDS_PAV_analysis.sh`](01_readsCDS_PAV_analysis.sh) - Whole-genome resequencing read mapping pipeline:

- Download SRA reads for multiple enset genotypes (BioProject [PRJNA344540](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA344540))
- Adapter trimming (Trim Galore)
- Read alignment to reference genomes (BWA-MEM2)
- Duplicate removal (Picard MarkDuplicates)
- Mapping quality assessment (Qualimap)
- Coverage calculation (BEDtools, deepTools)

### 3. Identify species-specific genes

[`02_identify_uniqueGenes.sh`](02_identify_uniqueGenes.sh) - Identification of genes unique to enset or *Musa* species based on PAV and CDS alignment results.
