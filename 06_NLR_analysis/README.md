# NLR Disease Resistance Gene Analysis

Identification and characterisation of nucleotide-binding leucine-rich repeat (NLR) disease resistance genes in *Ensete ventricosum* and *Musa* species genomes.

## Overview

NLR genes are a major class of intracellular immune receptors in plants. This module:

1. Identifies NLR loci directly in genome assemblies using motif-based detection (NLR-Annotator).
2. Predicts NLR proteins from gene models using HMM profiles and DIAMOND BLASTP.
3. Classifies NLRs into structural subclasses (CNL, TNL, NL, etc.).
4. Compares genomic NLR loci with predicted NLR gene models to detect non-annotated loci.
5. Reconstructs an NB-ARC domain phylogeny.

## Requirements

- `NLR-Annotator` (<https://github.com/steuernb/NLR-Annotator>)
- `MAST` / `FIMO` (MEME Suite)
- `HMMER` >= 3.3.2
- `DIAMOND` >= 2.0.15
- `InterProScan` >= 5.52
- `BLAST+` (tblastn, tblastx)
- `Clustal Omega` >= 1.2.1
- `trimAl` >= 1.4.1
- `RAxML` >= 8.2.12
- `Python` >= 3.9

Install conda-managed tools via the root [`environment.yml`](../environment.yml):

```bash
conda env create -f ../environment.yml
conda activate ev_mab_compgeno
```

NLR-Annotator and InterProScan require separate installation.

## Steps

### 1. Genomic NLR loci identification

[`NLR-annotator.sh`](NLR-annotator.sh) - Predict NLR loci directly from genome assemblies using [NLR-Annotator](https://github.com/steuernb/NLR-Annotator):

- Sequence chopping (ChopSequence)
- Motif parsing (NLR-Parser with MAST)
- NLR annotation and BED/GFF output

### 2. HMM-based NLR prediction from predicted proteins

[`NLR_predicted_prot.sh`](NLR_predicted_prot.sh) - Protein-level NLR prediction pipeline:

- Build custom HMM profiles from Pfam NLR domains (NB-ARC, TIR, RPW8, LRR) and the RefPlant NLR database
- DIAMOND BLASTP against curated NLR references
- HMMER search for NLR domain-containing proteins
- Fine-tuning with NB-ARC (PF00931) domain filter
- InterProScan functional annotation
- MAST/FIMO motif scanning

### 3. NLR class annotation

[`annotateClasses.py`](annotateClasses.py) - Classify predicted NLRs into classes (CNL, TNL, NL, CN, N, etc.) based on motif composition from NLR-Parser output.

Adapted from [Phillip Bayer's NLR summariser](https://gist.github.com/philippbayer/0052f5ad56121cd2252a1c5b90154ed1).

```bash
python annotateClasses.py --help
```

### 4. NLR loci vs NLR-encoding genes comparison

[`NLR_loci_vs_NLR_genes.sh`](NLR_loci_vs_NLR_genes.sh) - Map predicted NLR proteins back to genome assemblies (TBLASTN) and compare with NLR-Annotator genomic loci to identify non-overlapping NLR regions. Includes TBLASTX against NCBI nt for unmapped NLR loci.

### 5. Phylogenetic analysis

[`NLR_predicted_prot_phylogeny.sh`](NLR_predicted_prot_phylogeny.sh) - NB-ARC domain phylogeny:

- Multiple sequence alignment (Clustal Omega)
- Alignment trimming (trimAl, 90% gap threshold)
- Maximum likelihood tree (RAxML, PROTGAMMAAUTO, 1000 bootstraps)
