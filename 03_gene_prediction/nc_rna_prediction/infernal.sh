#!/bin/bash
# Infernal cmscan: whole-genome Rfam covariance model search (Rfam recommended workflow).
# Outputs per genome: *.cmscan.out, *.cmscan.tblout, *.cmscan.gff, *.rfam.fa
# Usage (SLURM array): sbatch --array=0-N infernal.sh  |  Interactive: bash infernal.sh <genome_prefix>

set -euo pipefail
trap 'rm -f "${TMPDIR:-/tmp}/$$".*' EXIT

# --- Helper functions ---
log()        { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err()        { log "ERROR: $*"; exit 1; }
check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool '$tool' not found in PATH"
    log "  $tool: $("$tool" -h 2>&1 | head -2 | tail -1 || true)"
}

# --- Configuration ---
MY_NUM_THREADS="${SLURM_CPUS_PER_TASK:-48}"
genotype="${SLURM_ARRAY_TASK_ID:-${1:?Usage: $0 <genome_prefix>}}"
DATA_DIR="${HOME}/data"
GENOME_FASTA="${DATA_DIR}/${genotype}.fna"
RFAM_CM="${HOME}/databases/Rfam/Rfam.cm"          # covariance model library
RFAM_CLANIN="${HOME}/databases/Rfam/Rfam.clanin"  # clan competition table
OUT_DIR="infernal"
VERSIONS_LOG="${OUT_DIR}/versions.log"
INFERNAL_BIN="${HOME}/apps/infernal-1.1.4/bin"

# --- Pre-flight checks ---
log "=== Infernal cmscan starting for genotype: ${genotype} ==="
for tool in "${INFERNAL_BIN}/cmscan" "${INFERNAL_BIN}/esl-sfetch" awk grep; do check_tool "$tool"; done
mkdir -p "$OUT_DIR"
[[ -f "$GENOME_FASTA" ]]  || err "Genome FASTA not found: ${GENOME_FASTA}"
[[ -f "$RFAM_CM" ]]       || err "Rfam covariance model not found: ${RFAM_CM}"
[[ -f "$RFAM_CLANIN" ]]   || err "Rfam clan table not found: ${RFAM_CLANIN}"
# Rfam.cm must be pre-pressed with 'cmpress Rfam.cm'; .i1m absence disables profile-HMM acceleration.
[[ -f "${RFAM_CM}.i1m" ]] \
    || err "Rfam.cm has not been pressed. Run: ${INFERNAL_BIN}/cmpress ${RFAM_CM}"

# Capture tool versions for reproducibility
{
    echo "=== Infernal ==="
    "${INFERNAL_BIN}/cmscan" -h 2>&1 | grep '# INFERNAL' | head -1 || true
    echo "=== Rfam.cm ==="
    head -3 "$RFAM_CM" || true
    echo "Run date:  $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Genotype:  ${genotype}"
    echo "Threads:   ${MY_NUM_THREADS}"
    echo "RFAM_CM:   ${RFAM_CM}"
    echo "RFAM_CLAN: ${RFAM_CLANIN}"
} >> "$VERSIONS_LOG"

# --- cmscan: search genome against all Rfam covariance models ---
# --nohmmonly: no HMM-only fallback (prevents false hits in low-complexity families)
# --rfam: Rfam heuristics for large genomes (required >100 Mb)
# --cut_ga: per-family gathering cutoffs set by Rfam curators
# --acc: report stable Rfam accessions (RF00001 etc.) instead of names
# --notextw: disable line-wrap for programmatic parsing
# --fmt 2 / --oclan / --oskip: extended tblout + clan competition (keeps top hit per locus)
# --clanin: clan membership table; --cpu: thread count
log "Running cmscan (Rfam, gathering cutoff, clan competition, ${MY_NUM_THREADS} CPUs)..."
"${INFERNAL_BIN}/cmscan" \
    --nohmmonly \
    --rfam \
    --cut_ga \
    --acc \
    --notextw \
    --fmt 2 \
    --oclan \
    --oskip \
    --clanin "$RFAM_CLANIN" \
    --cpu    "$MY_NUM_THREADS" \
    -o       "${OUT_DIR}/${genotype}.cmscan.out" \
    --tblout "${OUT_DIR}/${genotype}.cmscan.tblout" \
    "$RFAM_CM" \
    "$GENOME_FASTA"

N_HITS=$(grep -vc '^#' "${OUT_DIR}/${genotype}.cmscan.tblout" || true)
log "Rfam hits above gathering cutoff: ${N_HITS}"

# --- Index genome FASTA and extract hit subsequences ---
"${INFERNAL_BIN}/esl-sfetch" --index "$GENOME_FASTA"
# awk builds esl-sfetch key strings: "acc.seqname/start-end start end seqname"
grep -v '^#' "${OUT_DIR}/${genotype}.cmscan.tblout" \
    | awk '{ printf("%s.%s/%d-%d %d %d %s\n", $2, $4, $10, $11, $10, $11, $4); }' \
    | "${INFERNAL_BIN}/esl-sfetch" -Cf "$GENOME_FASTA" - \
    > "${OUT_DIR}/${genotype}.rfam.fa"
N_SEQS=$(grep -c '^>' "${OUT_DIR}/${genotype}.rfam.fa" || true)
log "Subsequences written to ${genotype}.rfam.fa: ${N_SEQS}"

# --- Convert tblout to GFF (cols: $4=seq $2=accession $10-$11=coords $17=score $12=strand) ---
grep -v '^#' "${OUT_DIR}/${genotype}.cmscan.tblout" \
    | awk 'BEGIN{OFS="\t"} {
        print $4, "infernal", $2, $10, $11, $17, $12, "."
    }' \
    > "${OUT_DIR}/${genotype}.cmscan.gff"

# --- Final summary ---
log "--- Summary: ${genotype} ---"
log "  Total cmscan hits (tblout lines): ${N_HITS}"
log "  Subsequences fetched (rfam.fa):   ${N_SEQS}"
log "  GFF records:  $(wc -l < "${OUT_DIR}/${genotype}.cmscan.gff")"
log "=== Infernal cmscan complete for genotype: ${genotype} ==="
