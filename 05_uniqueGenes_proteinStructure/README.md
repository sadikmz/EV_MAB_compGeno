# Protein Structure Prediction and Orthology Analysis of Species-Specific Genes

AlphaFold2-based protein structure prediction and OrthoFinder orthology analysis for genes identified as unique to *Ensete ventricosum* or *Musa* species in the comparative genomics step.

## Overview

1. **AlphaFold2 structure prediction** - Predict 3D protein structures for species-specific genes in monomer mode on a SLURM HPC cluster.
2. **Visualise AlphaFold2 results** - Generate per-protein coverage, predicted LDDT, and PAE plots.
3. **Orthology analysis** - Identify orthogroups across enset and *Musa* species using OrthoFinder with MSA-based phylogenetic inference.

## Requirements

- `AlphaFold2` (separate installation required; see <https://github.com/deepmind/alphafold>)
- `OrthoFinder` >= 2.5.5
- `MUSCLE` (for OrthoFinder MSA mode)
- `RAxML` >= 8.2.12
- `Python` >= 3.9 with `matplotlib`, `numpy`
- SLURM workload manager (for `run_alphafold.sbatch`)

Install conda-managed tools via the root [`environment.yml`](../environment.yml):

```bash
conda env create -f ../environment.yml
conda activate ev_mab_compgeno
```

## Steps

### 1. AlphaFold2 structure prediction

[`run_alphafold.sbatch`](run_alphafold.sbatch) - SLURM batch script to run AlphaFold2 (monomer mode) for protein structure prediction of species-specific genes.

Submit to a SLURM cluster with:

```bash
sbatch run_alphafold.sbatch
```

### 2. Visualise AlphaFold2 results

[`visualize_alphafold_results.py`](visualize_alphafold_results.py) - Python script (adapted from [VIBFold](https://github.com/jasperzuallaert/VIBFold)) to generate sequence coverage, predicted LDDT, and PAE plots from AlphaFold2 output.

[`visualize_alphafold_results.sh`](visualize_alphafold_results.sh) - Wrapper script to run the visualisation:

```bash
bash visualize_alphafold_results.sh
```

### 3. Orthology analysis

[`orthoFinderDir/run_orthofinder.sh`](orthoFinderDir/run_orthofinder.sh) - OrthoFinder analysis with MSA-based phylogenetic inference (MUSCLE + RAxML) to identify orthogroups across enset and *Musa* species.

```bash
bash orthoFinderDir/run_orthofinder.sh
```
