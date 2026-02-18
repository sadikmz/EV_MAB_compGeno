#!/bin/bash
# tRNAscan-SE: high-confidence eukaryotic tRNA prediction on a softmasked genome.
# Usage (interactive): bash tRNAscan-SE.sh <genome_prefix>
# Usage (SLURM array): sbatch --array=0-N tRNAscan-SE.sh
# Input:  <genome_prefix>.softmasked.fna   Outputs: tRNAscan/<genome_prefix>.*

set -euo pipefail
trap 'rm -f "${TMPDIR:-/tmp}/$$".*' EXIT

log()        { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err()        { log "ERROR: $*"; exit 1; }
check_tool() {
    command -v "$1" >/dev/null 2>&1 || err "Required tool '$1' not found in PATH"
    log "  $1: $("$1" --version 2>&1 | head -1)"
}

MY_NUM_THREADS="${SLURM_CPUS_PER_TASK:-28}"
export OMP_NUM_THREADS="$MY_NUM_THREADS"
genotype="${SLURM_ARRAY_TASK_ID:-${1:?Usage: $0 <genome_prefix>}}"
DATA_DIR="."
GENOME_FASTA="${DATA_DIR}/${genotype}.softmasked.fna"
OUT_DIR="tRNAscan"
VERSIONS_LOG="${OUT_DIR}/versions.log"

log "=== tRNAscan-SE starting: ${genotype} ==="
for tool in tRNAscan-SE awk; do check_tool "$tool"; done
mkdir -p "$OUT_DIR"
{ echo "=== tRNAscan-SE ==="; tRNAscan-SE --version 2>&1 || true
  echo "Run date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')  Genotype: ${genotype}  Threads: ${MY_NUM_THREADS}"
} >> "$VERSIONS_LOG"
[[ -f "$GENOME_FASTA" ]] || err "Genome FASTA not found: ${GENOME_FASTA}"

# -E: eukaryotic mode  --HQ: high-confidence (Infernal validated)
# --isotype: amino-acid identity annotation  --detail: intron/structure info
log "Running tRNAscan-SE (eukaryotic HQ mode, ${MY_NUM_THREADS} threads)..."
tRNAscan-SE -E --HQ --isotype --detail --threads "$MY_NUM_THREADS" \
    -o "${OUT_DIR}/${genotype}.tRNAscan.out" \
    -f "${OUT_DIR}/${genotype}.tRNAscan.ss" \
    -m "${OUT_DIR}/${genotype}.tRNAscan.stats" \
    -b "${OUT_DIR}/${genotype}.tRNAscan.bed" \
    "$GENOME_FASTA"

N_TRNA=$(awk 'NR>3 && NF>0' "${OUT_DIR}/${genotype}.tRNAscan.out" | grep -vc '^$' || true)
N_BED=$(wc -l < "${OUT_DIR}/${genotype}.tRNAscan.bed")
log "Summary: ${N_TRNA} tRNAs predicted  |  ${N_BED} BED intervals"
log "Isotype breakdown:"
awk 'NR>3 && NF>0 {print $5}' "${OUT_DIR}/${genotype}.tRNAscan.out" \
    | sort | uniq -c | sort -rn \
    | while read -r count iso; do log "  ${iso}: ${count}"; done

log "=== tRNAscan-SE complete: ${genotype} ==="
