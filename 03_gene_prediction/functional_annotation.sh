#!/bin/bash
set -euo pipefail

# Functional annotation of MAKER-predicted, TE-excluded proteins via:
#   1. BLASTp vs UniProt/Swiss-Prot  2. BLASTp vs UniProt/TrEMBL
#   3. BLASTp vs NCBI COG            4. BLASTp vs NCBI NR
#   5. EggNOG-mapper (synthesis / best functional annotation)
# All 4 BLAST searches run in parallel; EggNOG runs after they complete.

# --- Parameters (override via environment or edit here) --------------------
genotype="${GENOTYPE:-genotype_prefix}"
QUERY_FASTA="${QUERY_FASTA:-${genotype}.TE_excluded.fasta}"
UNIPROT_SPROT="${UNIPROT_SPROT:-${HOME}/data/uniprot/uniprot_sprot.fasta}"
UNIPROT_TREMBL="${UNIPROT_TREMBL:-${HOME}/data/uniprot/uniprot_trembl.fasta}"
COG_DB="${COG_DB:-${HOME}/data/cog_db/cog-20.fa}"
NR_DB="${NR_DB:-${HOME}/NCBI_NR/BLASTDB/nr}"
EGGNOG_DIR="${EGGNOG_DIR:-${HOME}/data/EggNOG}"
EGGNOG_MAPPER="${EGGNOG_MAPPER:-emapper.py}"
OUT_DIR="${OUT_DIR:-functional_annotation_out}"
cpus="${SLURM_CPUS_PER_TASK:-28}"
# Divide CPUs across 4 parallel BLAST jobs to avoid oversubscription
blast_cpus=$(( cpus / 4 )); (( blast_cpus < 1 )) && blast_cpus=1
EVALUE="1e-6"    # BLASTp e-value; EggNOG uses 1e-5 (diamond is more sensitive, so threshold can be relaxed)
BLAST_FMT=6

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

check_tool() {
    command -v "$1" >/dev/null 2>&1 || err "Required tool not found: $1"
    log "  $1: $("$1" --version 2>&1 | head -1)"
}

# Count annotated proteins in a BLAST tabular output (unique query hits)
summarise_blast() {
    local label="$1" outfile="$2"
    if [[ -f "${outfile}" ]]; then
        local n; n=$(cut -f1 "${outfile}" | sort -u | wc -l | tr -d ' ')
        log "  [${label}] Unique query proteins with hits: ${n}"
    else
        log "  [${label}] Output not found: ${outfile}"
    fi
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
log "=== Functional annotation pipeline ==="
log "Genotype  : ${genotype} | Query: ${QUERY_FASTA} | Output dir: ${OUT_DIR} | CPUs: ${cpus}"

for tool in blastp makeblastdb "${EGGNOG_MAPPER}"; do check_tool "${tool}"; done

[[ -f "${QUERY_FASTA}" ]] || err "Query FASTA not found: ${QUERY_FASTA}"
[[ -s "${QUERY_FASTA}" ]] || err "Query FASTA is empty: ${QUERY_FASTA}"
prot_count=$(grep -c '^>' "${QUERY_FASTA}" 2>/dev/null || echo 0)
log "Input proteins: ${prot_count}"

mkdir -p "${OUT_DIR}" "${OUT_DIR}/eggnog"

# ---------------------------------------------------------------------------
# Build BLAST databases (idempotent -- makeblastdb skips if already built)
# ---------------------------------------------------------------------------
log "--- Building BLAST databases ---"
[[ -f "${UNIPROT_SPROT}" ]]  || err "Swiss-Prot FASTA not found: ${UNIPROT_SPROT}"
[[ -f "${UNIPROT_TREMBL}" ]] || err "TrEMBL FASTA not found: ${UNIPROT_TREMBL}"
[[ -f "${COG_DB}" ]]         || err "COG FASTA not found: ${COG_DB}"

makeblastdb -in "${UNIPROT_SPROT}"  -dbtype prot -parse_seqids 2>&1 | tail -1 | { read l; log "  Swiss-Prot db: ${l}"; }
makeblastdb -in "${UNIPROT_TREMBL}" -dbtype prot -parse_seqids 2>&1 | tail -1 | { read l; log "  TrEMBL db: ${l}"; }
makeblastdb -in "${COG_DB}"         -dbtype prot -parse_seqids 2>&1 | tail -1 | { read l; log "  COG db: ${l}"; }
# NR is pre-built by NCBI; no makeblastdb step needed.

# ---------------------------------------------------------------------------
# Helper: single BLASTp search
# ---------------------------------------------------------------------------
run_blastp() {
    local label="$1" db="$2" outfile="$3"
    log "  [${label}] Starting BLASTp vs ${db}"
    blastp \
        -db              "${db}" \
        -query           "${QUERY_FASTA}" \
        -out             "${outfile}" \
        -evalue          "${EVALUE}" \
        -outfmt          "${BLAST_FMT}" \
        -max_target_seqs 1 \
        -seg             yes \
        -soft_masking    true \
        -lcase_masking \
        -max_hsps        1 \
        -num_threads     "${blast_cpus}"
    log "  [${label}] Finished -> ${outfile}"
}

# ---------------------------------------------------------------------------
# Run all 4 BLAST searches in parallel
# ---------------------------------------------------------------------------
log "--- Running BLAST searches (parallel) ---"
SPROT_OUT="${OUT_DIR}/${genotype}.vs_sprot.blastp.txt"
TREMBL_OUT="${OUT_DIR}/${genotype}.vs_trembl.blastp.txt"
COG_OUT="${OUT_DIR}/${genotype}.vs_cog.blastp.txt"
NR_OUT="${OUT_DIR}/${genotype}.vs_nr.blastp.txt"

run_blastp "Swiss-Prot" "${UNIPROT_SPROT}"  "${SPROT_OUT}"  &
run_blastp "TrEMBL"     "${UNIPROT_TREMBL}" "${TREMBL_OUT}" &
run_blastp "COG"        "${COG_DB}"         "${COG_OUT}"    &
run_blastp "NR"         "${NR_DB}"          "${NR_OUT}"     &

log "Waiting for all BLAST jobs to complete..."
wait
log "All BLAST searches complete."

# ---------------------------------------------------------------------------
# Summarise BLAST results
# ---------------------------------------------------------------------------
log "--- BLAST annotation summary ---"
summarise_blast "Swiss-Prot" "${SPROT_OUT}"
summarise_blast "TrEMBL"     "${TREMBL_OUT}"
summarise_blast "COG"        "${COG_OUT}"
summarise_blast "NR"         "${NR_OUT}"

# ---------------------------------------------------------------------------
# EggNOG-mapper -- synthesis functional annotation (diamond ultra-sensitive)
# ---------------------------------------------------------------------------
log "--- Running EggNOG-mapper ---"
DMND_DB="${EGGNOG_DIR}/eggnog_proteins.dmnd"
[[ -f "${DMND_DB}" ]] || err "EggNOG diamond db not found: ${DMND_DB}"
export EGGNOG_DATA_DIR="${EGGNOG_DIR}"

"${EGGNOG_MAPPER}" \
    --cpu                "${cpus}" \
    -m                   diamond \
    -o                   eggNOG \
    -i                   "${QUERY_FASTA}" \
    --itype              proteins \
    --sensmode           ultra-sensitive \
    --output_dir         "${OUT_DIR}/eggnog" \
    --dmnd_db            "${DMND_DB}" \
    --evalue             1e-5 \
    --override \
    --report_orthologs

log "EggNOG-mapper complete."

eggnog_annot="${OUT_DIR}/eggnog/eggNOG.emapper.annotations"
if [[ -f "${eggnog_annot}" ]]; then
    eg_hits=$(grep -c -v '^#' "${eggnog_annot}" 2>/dev/null || echo 0)
    log "  EggNOG annotated proteins: ${eg_hits}"
fi

log "=== Functional annotation complete. Outputs in: ${OUT_DIR} ==="
