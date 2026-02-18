# Gene Prediction and Functional Annotation

Protein-coding gene prediction, transcript assembly, and functional annotation for the *Ensete ventricosum* and *Musa* genome assemblies.

## Overview

The pipeline proceeds in the following steps:

1. **RNA-seq alignment** - Align raw RNA-seq reads to the repeat-masked genome (HISAT2).
2. **Transcript assembly** - Assemble transcripts (Trinity) and build splice-aware assemblies (PASA).
3. **Ab initio gene prediction** - Run BRAKER2 with GeneMark-ETP and AUGUSTUS using protein homology and RNA-seq evidence.
4. **Iterative gene prediction** - Train SNAP and run MAKER iteratively to integrate all evidence.
5. **Post-processing** - Filter gene models by AED score, assess completeness with BUSCO, and remove TE-associated models.
6. **Functional annotation** - Annotate predicted proteins.

## Requirements

- `HISAT2` >= 2.1.0
- `Trinity` >= 2.1.1
- `PASA` >= 2.4.1
- `BRAKER2` (with GeneMark-ETP and AUGUSTUS)
- `MAKER`
- `SNAP`
- `BUSCO` >= 5.2.1
- `TEsorter` (for TE-associated gene model filtering)
- `DIAMOND`, `eggNOG-mapper`, `InterProScan` (functional annotation)
- `R` >= 4.1 with `tidyverse` (for parsing and plotting)

Install conda-managed tools via the root [`environment.yml`](../environment.yml):

```bash
conda env create -f ../environment.yml
conda activate ev_mab_compgeno
```

BRAKER2, GeneMark-ETP, AUGUSTUS, MAKER, and PASA require separate installation following their respective documentation.

## Steps

### 1. Gene prediction

- [HISAT2 RNA-seq alignment](protein_coding_genes_prediction/hisat2.sh) - Align raw RNA-seq reads from landrace Mazia to the repeat-masked genome assembly using HISAT2 (v2.1.0).

- [Trinity transcript assembly](protein_coding_genes_prediction/trinity.sh) - Assemble RNA-seq transcripts *de novo* with Trinity (v2.1.1).

- [PASA pipeline](protein_coding_genes_prediction/pasa.sh) - Build splice-aware and unique transcript assemblies from Trinity output using the PASA pipeline (v2.4.1).

- [BRAKER2](protein_coding_genes_prediction/braker.sh) - Run BRAKER2 with GeneMark-ETP and AUGUSTUS to predict *ab initio* gene models supported by protein homology and/or RNA-seq evidence.

- [MAKER](protein_coding_genes_prediction/maker.sh) - Run MAKER iteratively to train SNAP gene models and integrate them with protein homology hints and assembled RNA-seq transcript evidence.

### 2. Extract and filter gene models, assess BUSCO scores

See [`collect_filter_maker_annotation.sh`](collect_filter_maker_annotation.sh).

Extracts predicted gene models from MAKER output, filters by Annotation Edit Distance (AED) score, and evaluates gene-space completeness with BUSCO.

### 3. Functional annotation of predicted proteins

See [`functional_annotation.sh`](functional_annotation.sh).

Performs functional annotation of filtered protein sequences using DIAMOND, eggNOG-mapper, and InterProScan.

### 4. Remove TE-associated gene models

Identifies and removes gene models that correspond to transposable element sequences:

- [`TEsorter_genemodels.sh`](TEsorter_genemodels.sh) - Run TEsorter on predicted proteins to flag TE-associated models.
- [`exclude_TE_inserted_gene_models.R`](exclude_TE_inserted_gene_models.R) - Parse TEsorter output and exclude TE-associated gene models from the final set.
- [`parse_plot_TEsorter_out.R`](parse_plot_TEsorter_out.R) - Plot TEsorter classification results.
