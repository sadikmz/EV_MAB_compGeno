#!/bin/bash
# hisat2.sh — RNA-seq alignment to softmasked genome for splice-site discovery
# Pipeline: (1) HISAT2 index, (2) build read lists, (3) align + sort BAM,
#           (4) validate BAM + mapping rate, (5) StringTie transcript assembly.
# Tools required: hisat2, samtools, stringtie

set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }
check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# --- Configuration — override via environment variables or edit here ---
genotype="${GENOTYPE:-genotype_prefix}"              # sample / accession identifier
softmasked_genome="${SOFTMASKED_GENOME:-softmasked.fna}"  # repeat-softmasked reference genome FASTA
RNASEQ_DIR="${RNASEQ_DIR:-${HOME}/data/rna_seq}"
RNASEQ_PREFIX="${RNASEQ_PREFIX:-rnaseq_lane_prefix}"    # set via env or edit here
MAX_INTRON_LEN="${MAX_INTRON_LEN:-160000}"           # max intron length (bp); monocot genomes ~100–200 kb
SUMMARY_FILE="${SUMMARY_FILE:-hisat2_alignment_summary.txt}"
MIN_MAPPING_RATE="${MIN_MAPPING_RATE:-50}"           # abort if overall mapping rate falls below this (%)
cpus="${SLURM_CPUS_PER_TASK:-28}"
export OMP_NUM_THREADS="${cpus}"
IFS=',' read -r -a LANES <<< "${LANES:-573,574,575}"

BASENAME="${softmasked_genome%.fna}"
INDEX_PREFIX="${softmasked_genome}.index"
OUTPUT_BAM="${BASENAME}.rnaseq.bam"
STRINGTIE_GFF="${BASENAME}_rnaseq.gff"
GENE_ABUND="gene_abundance.out"

# --- Tool checks ---
log "=== HISAT2 RNA-seq alignment pipeline ==="
for tool in hisat2 samtools stringtie; do check_tool "$tool"; done

# --- Input validation ---
log "Validating input files..."
[[ -f "$softmasked_genome" && -s "$softmasked_genome" ]] || err "Softmasked genome missing or empty: $softmasked_genome"
for lane in "${LANES[@]}"; do
    r1="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_1.fq.gz"
    r2="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_2.fq.gz"
    [[ -f "$r1" ]] || err "R1 reads not found for lane ${lane}: $r1"
    [[ -f "$r2" ]] || err "R2 reads not found for lane ${lane}: $r2"
done
log "All input files validated."

# --- Step 1: Build HISAT2 genome index (--large-index for genomes >4 Gb) ---
log "Step 1: Building HISAT2 genome index..."
if [[ -f "${INDEX_PREFIX}.1.ht2l" || -f "${INDEX_PREFIX}.1.ht2" ]]; then
    log "  Index already exists — skipping build."
else
    hisat2-build "$softmasked_genome" "$INDEX_PREFIX" -p "$cpus" --large-index
    log "  Index build complete."
fi
# Optional: extract known splice sites from a reference GFF to improve junction accuracy:
#   hisat2_extract_splice_sites.py reference.gff > known_splicesites.txt
#   then add --known-splicesite-infile known_splicesites.txt to the hisat2 call below.

# --- Step 2: Build comma-separated read lists from LANES array ---
log "Step 2: Preparing read lists for ${#LANES[@]} lanes..."
R1_LIST=""; R2_LIST=""
for lane in "${LANES[@]}"; do
    r1="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_1.fq.gz"
    r2="${RNASEQ_DIR}/${RNASEQ_PREFIX}_${lane}_2.fq.gz"
    R1_LIST="${R1_LIST:+${R1_LIST},}${r1}"
    R2_LIST="${R2_LIST:+${R2_LIST},}${r2}"
done
log "  R1 files: $R1_LIST"
log "  R2 files: $R2_LIST"

# --- Step 3: Align RNA-seq reads ---
# --dta: optimised for StringTie (do NOT combine with --dta-cufflinks)
# --no-mixed / --no-discordant: keep only properly paired alignments
# --very-sensitive: highest sensitivity preset; --max-intronlen: plant introns rarely >160 kb
log "Step 3: Aligning RNA-seq reads with HISAT2..."
hisat2 \
    -q \
    -x "$INDEX_PREFIX" \
    -p "$cpus" \
    -1 "$R1_LIST" \
    -2 "$R2_LIST" \
    --dta \
    --max-intronlen "$MAX_INTRON_LEN" \
    --no-mixed \
    --very-sensitive \
    --no-discordant \
    --summary-file "$SUMMARY_FILE" \
    | samtools view -@ "$cpus" -bS \
    | samtools sort -@ "$cpus" -o "$OUTPUT_BAM" -T hisat2_tmp
log "  Alignment and sort complete: $OUTPUT_BAM"

# --- Step 4: Validate BAM and mapping rate ---
log "Step 4: Validating BAM file and mapping rate..."
[[ -f "$OUTPUT_BAM" && -s "$OUTPUT_BAM" ]] || err "BAM missing or empty: $OUTPUT_BAM"
samtools index "$OUTPUT_BAM"
FLAGSTAT_OUT="${BASENAME}.rnaseq.flagstat.txt"
samtools flagstat "$OUTPUT_BAM" | tee "$FLAGSTAT_OUT"
log "  Flagstat report written to $FLAGSTAT_OUT"

# Parse overall mapping rate from HISAT2 summary ("XX.XX% overall alignment rate")
if [[ -f "$SUMMARY_FILE" ]]; then
    MAP_RATE=$(grep "overall alignment rate" "$SUMMARY_FILE" \
        | grep -oP '[0-9]+\.[0-9]+' | head -1 || echo "0")
    MAP_INT=${MAP_RATE%.*}
    log "  Overall mapping rate: ${MAP_RATE}%"
    if (( MAP_INT < MIN_MAPPING_RATE )); then
        err "Mapping rate ${MAP_RATE}% is below the minimum threshold of ${MIN_MAPPING_RATE}%. Check read quality or genome assembly."
    fi
    log "  Mapping rate OK (>= ${MIN_MAPPING_RATE}%)."
fi

# --- Step 5: Transcript assembly with StringTie ---
# Produces GTF/GFF of assembled transcripts + per-gene abundance for PASA/MAKER evidence
log "Step 5: Assembling transcripts with StringTie..."
stringtie "$OUTPUT_BAM" -o "$STRINGTIE_GFF" -p "$cpus" -A "$GENE_ABUND"
[[ -f "$STRINGTIE_GFF" ]] || err "StringTie GFF not produced: $STRINGTIE_GFF"
TRANSCRIPT_COUNT=$(grep -c $'\ttranscript\t' "$STRINGTIE_GFF" || true)
log "  StringTie assembled ${TRANSCRIPT_COUNT} transcripts -> $STRINGTIE_GFF"

log "=== HISAT2 pipeline complete ==="
log "  BAM:           $OUTPUT_BAM"
log "  BAM index:     ${OUTPUT_BAM}.bai"
log "  Flagstat:      $FLAGSTAT_OUT"
log "  Alignment log: $SUMMARY_FILE"
log "  StringTie GFF: $STRINGTIE_GFF"
log "  Gene abundance: $GENE_ABUND"
