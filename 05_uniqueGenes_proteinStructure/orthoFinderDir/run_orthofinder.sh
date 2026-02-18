#!/bin/bash
# run_orthofinder.sh — OrthoFinder orthogroup inference across EV/MAB and outgroup proteomes.
# Steps: DIAMOND all-vs-all → MCL clustering → MUSCLE MSA → RAxML gene trees → orthologue assignment.

set -euo pipefail

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }

check_tool() {
    local tool="$1"
    if ! command -v "${tool}" &>/dev/null; then
        err "Required tool not found on PATH: ${tool}"
        err "Ensure the relevant module is loaded and the tool is installed."
    fi
}

# ---------------------------------------------------------------------------
# Configurable parameters
# ---------------------------------------------------------------------------
PROT_DIR="${PROT_DIR:-/path/to/ev_mab_prots}"
INPUT_DIR="${INPUT_DIR:-${PROT_DIR}}"
PATH_GENOMES="${PATH_GENOMES:-/path/to/genome_dir}"
QPROT_DIR="${QPROT_DIR:-/path/to/EV_MAB_novel_genes_db}"
DMND_BLASTDB="${DMND_BLASTDB:-/path/to/NCBI_NR}"

# -t: DIAMOND all-vs-all threads (I/O-heavy); -a: MCL/tree threads (~1/10 of -t, floor 1)
threads="${SLURM_NTASKS_PER_NODE:-80}"
algo_threads=$(( threads / 10 ))
algo_threads=$(( algo_threads < 1 ? 1 : algo_threads ))

# ---------------------------------------------------------------------------
# Module environment
# ---------------------------------------------------------------------------
module purge
# RAxML hybrid (MPI+PThreads, avx2); PROTGAMMAAUTO selects best model per gene tree.
module load RAxML/8.2.12-gompi-2020a-hybrid-avx2

# ---------------------------------------------------------------------------
# Tool checks and version logging
# ---------------------------------------------------------------------------
for tool in orthofinder muscle raxmlHPC-PTHREADS-AVX2; do check_tool "${tool}"; done

log "=== Tool versions ==="
log "OrthoFinder: $(orthofinder --version 2>&1 | head -1)"
log "RAxML:       $(raxmlHPC-PTHREADS-AVX2 -v 2>&1 | grep -m1 'RAxML')"
log "MUSCLE:      $(muscle --version 2>&1 | head -1 || muscle -version 2>&1 | head -1)"

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
log "=== Input validation ==="
[[ -d "${INPUT_DIR}" ]] || { err "Input proteome directory not found: ${INPUT_DIR}"; err "Set PROT_DIR or INPUT_DIR to the directory containing per-species protein FASTAs."; exit 1; }

n_species=$(find "${INPUT_DIR}" -maxdepth 1 \( -name "*.fa" -o -name "*.faa" \) | wc -l)
[[ "${n_species}" -gt 0 ]] || { err "No .fa or .faa files found in ${INPUT_DIR}"; err "OrthoFinder requires one protein FASTA file per species."; exit 1; }

log "Found ${n_species} species proteome file(s) in ${INPUT_DIR}"
log "Genome directory (provenance):    ${PATH_GENOMES}"
log "Novel/query proteome directory:   ${QPROT_DIR}"
log "DIAMOND NR database (provenance): ${DMND_BLASTDB}"

# ---------------------------------------------------------------------------
# Run OrthoFinder
# ---------------------------------------------------------------------------
log "=== Starting OrthoFinder ==="
log "Input: ${INPUT_DIR} | -t ${threads} | -a ${algo_threads} | -M msa | -A muscle | -T raxml"

orthofinder \
    -f "${INPUT_DIR}" \
    -t "${threads}" \
    -a "${algo_threads}" \
    -M msa \
    -A muscle \
    -T raxml

# ---------------------------------------------------------------------------
# Post-run summary
# ---------------------------------------------------------------------------
log "=== Post-run summary ==="

# OrthoFinder writes results to a timestamped subdirectory inside INPUT_DIR.
results_dir=$(find "${INPUT_DIR}" -maxdepth 2 -type d -name "OrthoFinder" | sort | tail -1)
[[ -n "${results_dir}" ]] || { err "Could not locate OrthoFinder results directory under ${INPUT_DIR}"; exit 1; }

results_run=$(find "${results_dir}" -maxdepth 1 -type d -name "Results_*" | sort | tail -1)
log "Results directory: ${results_run}"

orthogroups_tsv="${results_run}/Orthogroups/Orthogroups.tsv"
if [[ -f "${orthogroups_tsv}" ]]; then
    n_orthogroups=$(( $(wc -l < "${orthogroups_tsv}") - 1 ))
    log "Orthogroups inferred: ${n_orthogroups}"
else
    log "Orthogroups.tsv not found; OrthoFinder may not have completed successfully."
fi

stats_tsv="${results_run}/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"
if [[ -f "${stats_tsv}" ]]; then
    log "Per-species statistics:"
    column -t "${stats_tsv}" | head -20 | while IFS= read -r line; do log "  ${line}"; done
fi

log "=== OrthoFinder run complete ==="
