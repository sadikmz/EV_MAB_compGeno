#!/bin/bash
# trinity.sh — Quality trimming (Trim Galore) + de novo assembly (Trinity >= v2.15.1)
# Pipeline: (1) trim per lane, (2) validate trimmed files, (3) Trinity assembly,
#           (4) validate output + TrinityStats.pl metrics.

set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }
check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# --- Configuration — override via environment variables or edit here ---
TRINITY_HOME="${TRINITY_HOME:-${HOME}/apps/trinityrnaseq-v2.15.1}"
RNASEQ_DIR="${RNASEQ_DIR:?RNASEQ_DIR must be set to the directory containing raw read files}"
IFS=',' read -r -a LANES <<< "${LANES:-573,574,575}"
RNASEQ_PREFIX="${RNASEQ_PREFIX:-rnaseq_lane_prefix}"    # set via env or edit here
SS_LIB_TYPE="${SS_LIB_TYPE:-RF}"        # RF = dUTP / reverse-stranded libraries
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-100}" # discard contigs shorter than this (bp)
MAX_MEMORY="${MAX_MEMORY:-100G}"        # RAM for Jellyfish / Inchworm stages
TRIM_QUALITY="${TRIM_QUALITY:-30}"      # Phred quality threshold for Trim Galore
OUTPUT_DIR="${OUTPUT_DIR:-trinity_out}"
cpus="${SLURM_CPUS_PER_TASK:-28}"
TRIM_DIR="${TRIM_DIR:-trim_out}"
TRINITY_FASTA="${OUTPUT_DIR}/Trinity.fasta"

# --- Tool checks ---
log "=== Trinity de novo transcriptome assembly pipeline ==="
check_tool trim_galore
[[ -x "${TRINITY_HOME}/Trinity" ]] || err "Trinity not found at: ${TRINITY_HOME}/Trinity"
log "  Trinity: $("${TRINITY_HOME}/Trinity" --version 2>&1 | head -1)"
TRINITY_STATS="${TRINITY_HOME}/util/TrinityStats.pl"
[[ -x "$TRINITY_STATS" ]] || err "TrinityStats.pl not found at: $TRINITY_STATS"

# --- Validate raw reads ---
log "Validating raw read files..."
for lane in "${LANES[@]}"; do
    r1="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_1.fq.gz"
    r2="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_2.fq.gz"
    [[ -f "$r1" && -s "$r1" ]] || err "R1 reads missing or empty for lane ${lane}: $r1"
    [[ -f "$r2" && -s "$r2" ]] || err "R2 reads missing or empty for lane ${lane}: $r2"
done

# --- Step 1: Trim Galore (per lane) ---
# --quality 30: trim 3' bases below Phred 30; --paired: both mates must pass length filter
log "Step 1: Quality trimming with Trim Galore (Q${TRIM_QUALITY})..."
mkdir -p "$TRIM_DIR"
for lane in "${LANES[@]}"; do
    r1="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_1.fq.gz"
    r2="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_2.fq.gz"
    log "  Trimming lane ${lane}..."
    trim_galore --cores "$cpus" --quality "$TRIM_QUALITY" --output_dir "$TRIM_DIR" --paired "$r1" "$r2"
done

# --- Step 2: Validate trimmed files (val_1/val_2 naming from Trim Galore) ---
log "Step 2: Validating trimmed read files..."
LEFT_READS=(); RIGHT_READS=()
for lane in "${LANES[@]}"; do
    val1="${TRIM_DIR}/${RNASEQ_PREFIX}_${lane}_1_val_1.fq.gz"
    val2="${TRIM_DIR}/${RNASEQ_PREFIX}_${lane}_2_val_2.fq.gz"
    [[ -f "$val1" && -s "$val1" ]] || err "Trimmed R1 missing or empty for lane ${lane}: $val1"
    [[ -f "$val2" && -s "$val2" ]] || err "Trimmed R2 missing or empty for lane ${lane}: $val2"
    LEFT_READS+=("$val1"); RIGHT_READS+=("$val2")
done
LEFT_LIST=$(IFS=,; echo "${LEFT_READS[*]}")
RIGHT_LIST=$(IFS=,; echo "${RIGHT_READS[*]}")

# --- Step 3: Trinity assembly ---
# --normalize_reads: caps ~50x coverage, reducing memory while preserving diversity
# --SS_lib_type RF: dUTP libraries (R=first read antisense, F=second read sense)
# --full_cleanup: remove large intermediates; --monitor: print resource usage every 60s
log "Step 3: Running Trinity assembly..."
export PERL5LIB=""   # prevent system Perl modules conflicting with Trinity
"${TRINITY_HOME}/Trinity" \
    --seqType fq \
    --left  "$LEFT_LIST" \
    --right "$RIGHT_LIST" \
    --SS_lib_type "$SS_LIB_TYPE" \
    --min_contig_length "$MIN_CONTIG_LEN" \
    --max_memory "$MAX_MEMORY" \
    --CPU "$cpus" \
    --normalize_reads \
    --full_cleanup \
    --monitor \
    --output "$OUTPUT_DIR"

# --- Step 4: Validate output and report statistics ---
log "Step 4: Validating Trinity output and computing assembly statistics..."
[[ -f "$TRINITY_FASTA" ]] || err "Trinity FASTA not produced: $TRINITY_FASTA"
[[ -s "$TRINITY_FASTA" ]] || err "Trinity FASTA is empty: $TRINITY_FASTA"
CONTIG_COUNT=$(grep -c '^>' "$TRINITY_FASTA" || true)
log "  Total assembled contigs: $CONTIG_COUNT"
log "  Assembly statistics (TrinityStats.pl):"
"$TRINITY_STATS" "$TRINITY_FASTA" 2>&1 | while IFS= read -r line; do log "    $line"; done

log "=== Trinity pipeline complete ==="
log "  Assembly FASTA: $TRINITY_FASTA"
log "  Total contigs:  $CONTIG_COUNT"
log "  Trimmed reads:  $TRIM_DIR/"
