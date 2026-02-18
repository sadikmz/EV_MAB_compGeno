# Comparative Genomics of *Ensete ventricosum* and *Musa* Species

Scripts and pipelines accompanying the manuscript:

> **What makes a banana false? How the genome of Ethiopian orphan staple *Ensete ventricosum* differs from the banana A and B sub-genomes**
>
> Muzemil S, Paul P, Baxter L, Dominguez-Ferreras A, Sahu SK, Van Deynze A, Mai G, Yemataw Z, Tesfaye K, Ntoukakis V, Studholme DJ, Grant M (2026).

## Overview

This repository contains the bioinformatics pipelines used for *de novo* genome assembly quality assessment, repeat identification and annotation, gene prediction, comparative genomic analysis, protein structure prediction, and NLR disease resistance gene characterisation of *Ensete ventricosum* (enset) landraces Mazia and Bedadeti, compared with banana progenitors *Musa acuminata* (A genome) and *Musa balbisiana* (B genome).

## Getting Started

### Requirements

All dependencies are specified in `environment.yml`. To create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate ev_mab_compgeno
```

Some tools (e.g. AlphaFold2, BRAKER2, NLR-Annotator, InterProScan) require separate installation; follow the links in the relevant subdirectory READMEs.

## Repository Structure

| Directory | Description |
|-----------|-------------|
| [`01_genome_qc_and_profiling`](01_genome_qc_and_profiling/) | Genome size estimation (jellyfish, GenomeScope2, FindGSE), ploidy assessment (Smudgeplot), and assembly QC (gfastats, BUSCO) |
| [`02_repeat_identification`](02_repeat_identification/) | Structural repeat identification (RepeatMasker, RepeatModeler, MITE-Hunter, EDTA, PILER, TRF) and annotation (TEsorter) |
| [`03_gene_prediction`](03_gene_prediction/) | Gene prediction (BRAKER2, MAKER), RNA-seq alignment (HISAT2), transcript assembly (Trinity, PASA), and functional annotation |
| [`04_comparative_genomics`](04_comparative_genomics/) | Presence/absence variation (PAV) analysis via read mapping and CDS alignment; identification of species-specific genes |
| [`05_uniqueGenes_proteinStructure`](05_uniqueGenes_proteinStructure/) | AlphaFold2 protein structure prediction for species-specific genes; OrthoFinder orthology analysis |
| [`06_NLR_analysis`](06_NLR_analysis/) | NLR disease resistance gene identification (NLR-Annotator, HMM-based prediction), class annotation, phylogenetic analysis |

## Key Tools and Dependencies

- **Genome QC**: jellyfish, GenomeScope2, FindGSE, Smudgeplot (FastK, PloidyPlot), gfastats, BUSCO
- **Repeat analysis**: RepeatMasker, RepeatModeler, MITE-Hunter, TRF, PILER, CD-HIT, EDTA, TEsorter
- **Gene prediction**: BRAKER2 (GeneMark-ETP, AUGUSTUS), MAKER, SNAP, HISAT2, Trinity, PASA
- **Comparative genomics**: BWA-MEM2, SAMtools, BEDtools, Picard, Qualimap, deepTools
- **Protein structure**: AlphaFold2, OrthoFinder
- **NLR analysis**: NLR-Annotator, HMMER, InterProScan, DIAMOND, MAST/FIMO
- **Phylogenetics**: Clustal Omega, trimAl, RAxML

## Data Availability

Genome assemblies and raw sequencing data are available under NCBI BioProject [PRJNA838108](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA838108).

## Citation

If you use these scripts, please cite the accompanying manuscript (details above).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
