#!/bin/bash
# run_alphafold.sbatch — SLURM batch script to run AlphaFold2 on a single protein query.
# Based on: https://elearning.vib.be/courses/alphafold/lessons/alphafold-on-the-hpc/topic/setting-up-the-job-script/
# Usage: sbatch run_alphafold.sbatch
# Override defaults via environment: PROTEIN=myprotein.faa PRESET=multimer sbatch run_alphafold.sbatch

set -euo pipefail

echo "hostname=$(hostname)"
echo "pwd=$(pwd)"
echo "USER=${USER}"
date +"%d-%b-%y %T"
# ====================

module load AlphaFold

# AlphaFold database directory — override with ALPHAFOLD_DATA_DIR env var if needed.
ALPHAFOLD_DATA_DIR="${ALPHAFOLD_DATA_DIR:-/home/data/bioinf_resources/Alphafold/20231026}"
export ALPHAFOLD_DATA_DIR

# Configuration — override via environment variables before submitting.
PROTEIN="${PROTEIN:-rbd.faa}"
PRESET="${PRESET:-monomer_ptm}"    # monomer_ptm or multimer
NAME="$(basename "${PROTEIN}" .faa)"
OUTDIR="./${NAME}"

[[ -f "${PROTEIN}" ]] || { echo "ERROR: FASTA not found: ${PROTEIN}" >&2; exit 1; }
mkdir -p "${OUTDIR}"

alphafold \
    --fasta_paths="${PROTEIN}" \
    --max_template_date=2999-01-01 \
    --db_preset=full_dbs \
    --output_dir="${OUTDIR}" \
    --model_preset="${PRESET}"

echo "Done SLURM job: ${SLURM_JOBID:-<interactive>}"
echo "Done running $0"
date +"%d-%b-%y %T"
