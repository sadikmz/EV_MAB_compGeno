#!/bin/bash
# rRNA homology search: BLASTN against curated plant rRNA templates
# (Arabidopsis 5S/5.8S/18S; Oryza sativa 28S). Results converted to sorted,
# merged BED intervals for use in the cascaded ncRNA pipeline (ncRNA_summary.sh).
# Usage: sbatch rRNA_homology_search.sh  OR  bash rRNA_homology_search.sh

set -euo pipefail
trap 'rm -f "${TMPDIR:-/tmp}/$$".*' EXIT

log()        { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err()        { log "ERROR: $*"; exit 1; }
check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool '$tool' not found in PATH"
    log "  $tool: $("$tool" -version 2>&1 | head -1 || "$tool" --version 2>&1 | head -1 || true)"
}

# --- Configuration — edit paths / parameters as needed ---
MY_NUM_THREADS="${SLURM_CPUS_PER_TASK:-28}"
DATA_DIR="${HOME}/data"                  # directory containing *.fna genome assemblies
rRNA_DB="${HOME}/data/rRNA/rRNA.fna"    # FASTA of curated rRNA template sequences
OUT_DIR="homology_search"
VERSIONS_LOG="${OUT_DIR}/versions.log"
EVALUE="1e-5"       # permissive e-value; rRNA fragments can be short/partially degraded
PERC_IDENTITY="90"  # ≥90% identity (plant 18S/28S ~95–99%; 5S can diverge at flanks)
WORD_SIZE="11"      # blastn default; more sensitive than megablast for diverged fragments
TASK="blastn"       # more sensitive than megablast for rRNA fragment detection

# --- Pre-flight checks ---
log "=== rRNA homology search (BLASTN) starting ==="
for tool in blastn makeblastdb bedtools awk; do check_tool "$tool"; done
mkdir -p "$OUT_DIR"
[[ -f "$rRNA_DB" ]] || err "rRNA reference FASTA not found: ${rRNA_DB}"
[[ -d "$DATA_DIR" ]] || err "Data directory not found: ${DATA_DIR}"

# Log tool versions and run parameters for reproducibility
{
    echo "=== BLAST+ ==="; blastn -version 2>&1 | head -2 || true
    echo "=== bedtools ==="; bedtools --version 2>&1 | head -1 || true
    echo "Run date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "rRNA_DB: ${rRNA_DB}  EVALUE: ${EVALUE}  PERC_ID: ${PERC_IDENTITY}  TASK: ${TASK}"
} >> "$VERSIONS_LOG"

# --- Build BLAST database (skip if .nhr index already exists) ---
if [[ ! -f "${rRNA_DB}.nhr" ]]; then
    log "Building BLAST nucleotide database from: ${rRNA_DB}"
    makeblastdb -in "$rRNA_DB" -input_type fasta -dbtype nucl -logfile "${OUT_DIR}/makeblastdb.log"
    log "Database build complete."
else
    log "BLAST database index found (${rRNA_DB}.nhr); skipping makeblastdb."
fi

# --- Search each genome assembly ---
shopt -s nullglob
GENOMES=("${DATA_DIR}"/*.fna)
[[ ${#GENOMES[@]} -gt 0 ]] || err "No *.fna files found in ${DATA_DIR}"

for genome_path in "${GENOMES[@]}"; do
    BASENAME=$(basename "$genome_path" .fna)
    log "Processing: ${BASENAME}"

    # BLASTN: -outfmt 6 tabular; -num_alignments 1 keeps best hit per query region
    blastn \
        -task         "$TASK" \
        -db           "$rRNA_DB" \
        -query        "$genome_path" \
        -outfmt       6 \
        -evalue       "$EVALUE" \
        -perc_identity "$PERC_IDENTITY" \
        -word_size    "$WORD_SIZE" \
        -num_threads  "$MY_NUM_THREADS" \
        -num_alignments 1 \
        -out          "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.out"

    # Convert fmt6 to 0-based BED; swap coords for minus-strand hits (qstart > qend)
    awk 'BEGIN{OFS="\t"} {
        qs = $7 - 1; qe = $8
        if (qs > qe) { tmp = qs; qs = qe; qe = tmp }
        print $1, qs, qe, $2
    }' "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.out" \
        > "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.bed"

    # Sort, normalise columns, then merge overlapping rRNA intervals
    bedtools sort -i "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.bed" \
        > "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.bed"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
        "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.bed" \
        > "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.1.bed"
    # -c 4 -o distinct: collapse overlapping hit labels to preserve rRNA subtype info
    bedtools merge \
        -i "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.1.bed" \
        -c 4 -o distinct \
        > "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.merged.bed"

    N_RAW=$(wc -l < "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.out")
    N_MERGED=$(wc -l < "${OUT_DIR}/${BASENAME}.rRNA.e05.blastn.sorted.merged.bed")
    log "  ${BASENAME}: ${N_RAW} raw BLASTN hits → ${N_MERGED} merged rRNA loci"
done

log "=== rRNA homology search complete ==="
