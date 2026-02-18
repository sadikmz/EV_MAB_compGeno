#!/bin/bash
# TEsorter: classify TE-related domains in MAKER-predicted gene model proteins.
# Input: *.faa protein files (EV_mazia, EV_bedadeti, Musa_acuminata,
#        Musa_balbisiana, Ensete_glaucum); EG/Musa from banana-genome-hub.southgreen.fr

set -euo pipefail

INPUT_DIR="${INPUT_DIR:-${HOME}/data}"
OUTPUT_DIR="${OUTPUT_DIR:-tesorter_out}"
HMM_DB="${HMM_DB:-rexdb-plant}"
PASS2_RULE="${PASS2_RULE:-80-80-80}"
cpus="${SLURM_NTASKS_PER_NODE:-48}"
# Optional: GENOME_GFF3="${GENOME_GFF3:-}"  (restrict to gene model loci)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }
check_tool() {
    command -v "$1" >/dev/null 2>&1 || err "Required tool not found: $1"
    log "  $1: $("$1" --version 2>&1 | head -1)"
}

log "=== TEsorter gene model classification ==="
log "Input: ${INPUT_DIR}  Output: ${OUTPUT_DIR}  DB: ${HMM_DB}  CPUs: ${cpus}"

check_tool TEsorter

shopt -s nullglob; input_files=("${INPUT_DIR}"/*.faa); shopt -u nullglob
[[ ${#input_files[@]} -gt 0 ]] || err "No .faa files found in ${INPUT_DIR}"
log "Found ${#input_files[@]} input protein file(s)"
mkdir -p "${OUTPUT_DIR}"

for faa_file in "${input_files[@]}"; do
    BASENAME=$(basename "${faa_file}" .faa)
    PREFIX="${OUTPUT_DIR}/${BASENAME}_TEsorter.${HMM_DB}"
    log "--- Processing: ${BASENAME} ---"
    if [[ ! -s "${faa_file}" ]]; then log "WARNING: ${faa_file} is empty -- skipping"; continue; fi
    log "  Input proteins: $(grep -c '^>' "${faa_file}" 2>/dev/null || echo 0)"

    TEsorter "${faa_file}" \
        --hmm        "${HMM_DB}" \
        --seq-type   prot \
        --processors "${cpus}" \
        --pass2-rule "${PASS2_RULE}" \
        --prefix     "${PREFIX}"
    # Optional: add --genome-anno "${GENOME_GFF3}" to restrict to gene model loci

    for tsv in cls dom; do
        f="${PREFIX}.${tsv}.tsv"
        [[ -f "$f" ]] \
            && log "  ${tsv}.tsv hits: $(grep -c '.' "$f" 2>/dev/null || echo 0)" \
            || log "  WARNING: ${tsv}.tsv not produced for ${BASENAME}"
    done
    log "  Done: ${BASENAME}"
done

log "=== TEsorter complete. Outputs in: ${OUTPUT_DIR} ==="
