#!/bin/bash
# pasa.sh — Transcript alignment and assembly using PASA (>= v2.5.1).
# Pipeline: seqclean → BLAT+GMAP alignment → SQLite PASA assembly database.
# Ref: https://github.com/PASApipeline/PASApipeline

set -euo pipefail

# --- Logging helpers ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

# --- Configuration — override via environment variables or edit here ---
genotype="${GENOTYPE:-genotype_prefix}"
genome="${GENOME:-${genotype}.fna}"
trinity_raw="${TRINITY_RAW:-trinity_out/Trinity.fasta}"
trinity_cleaned="${TRINITY_CLEANED:-data/Trinity.cleaned.fasta}"
trans_gtf="${TRANS_GTF:-${genotype}_rnaseq.gff}"
PASA_DIR="${PASA_DIR:-${HOME}/apps/PASApipeline-v2.5.1}"
PASA_CONFIG="${PASA_CONFIG:-${PASA_DIR}/pasa_conf/pasa.config}"
DATABASE="${DATABASE:-${PWD}/pasa_${genotype}/pasa.sqlite}"  # persistent, deterministic; survives restarts
SEQCLEAN_CPUS="${SEQCLEAN_CPUS:-48}"
SEQCLEAN_MAX_N="${SEQCLEAN_MAX_N:-10000}"
cpus="${SLURM_CPUS_PER_TASK:-28}"

# --- Tool checks ---
log "=== PASA transcript alignment and assembly pipeline ==="
SEQCLEAN_BIN="${PASA_DIR}/bin/seqclean"
PASA_BIN="${PASA_DIR}/Launch_PASA_pipeline.pl"
[[ -x "$SEQCLEAN_BIN" ]] || err "seqclean not found at: $SEQCLEAN_BIN"
log "  seqclean: $("$SEQCLEAN_BIN" 2>&1 | head -1 || true)"
[[ -x "$PASA_BIN" ]]     || err "Launch_PASA_pipeline.pl not found at: $PASA_BIN"
log "  PASA: $("$PASA_BIN" --version 2>&1 | head -1 || true)"
[[ -f "$PASA_CONFIG" ]]  || err "PASA config not found: $PASA_CONFIG"

# --- Input validation ---
log "Validating input files..."
[[ -f "$trinity_raw"  && -s "$trinity_raw"  ]] || err "Trinity FASTA not found or empty: $trinity_raw"
[[ -f "$genome"       && -s "$genome"       ]] || err "Genome FASTA not found or empty: $genome"
[[ -f "$trans_gtf"    && -s "$trans_gtf"    ]] || err "StringTie GFF not found or empty: $trans_gtf"
RAW_COUNT=$(grep -c '^>' "$trinity_raw" || true)
log "  Input Trinity transcripts: $RAW_COUNT — all inputs validated."

# --- Step 1: seqclean (removes low-complexity, excess-N, and vector sequences) ---
log "Step 1: Cleaning Trinity transcripts with seqclean..."
if [[ -f "$trinity_cleaned" && -s "$trinity_cleaned" ]]; then
    log "  Cleaned FASTA already exists — skipping seqclean: $trinity_cleaned"
else
    "$SEQCLEAN_BIN" "$trinity_raw" -c "$SEQCLEAN_CPUS" -o "$trinity_cleaned" -n "$SEQCLEAN_MAX_N"
fi
[[ -f "$trinity_cleaned" ]] || err "seqclean did not produce output: $trinity_cleaned"
[[ -s "$trinity_cleaned" ]] || err "seqclean output is empty — all transcripts were removed: $trinity_cleaned"
CLEANED_COUNT=$(grep -c '^>' "$trinity_cleaned" || true)
log "  Transcripts after cleaning: ${CLEANED_COUNT} (removed $((RAW_COUNT - CLEANED_COUNT)))"

# --- Step 2: PASA alignment and assembly (BLAT + GMAP, SQLite DB, TransDecoder) ---
# -C: create new DB  -R: run alignment/assembly  --TRANSDECODER: ORF calling within PASA
# --stringent_alignment_overlap 0.30: require >=30% overlap to prevent spurious merging
log "Step 2: Running PASA alignment and assembly..."
# Patch DATABASE into a run-specific config copy (--PASACONF is not a valid PASA flag).
RUNTIME_CONFIG="pasa_runtime_${genotype}.config"
sed "s|^DATABASE=.*|DATABASE=${DATABASE}|" "$PASA_CONFIG" > "$RUNTIME_CONFIG"
mkdir -p "$(dirname "$DATABASE")"
"$PASA_BIN" \
    -c "$RUNTIME_CONFIG" \
    -C -R \
    -g "$genome" \
    -t "$trinity_cleaned" \
    --trans_gtf "$trans_gtf" \
    --TRANSDECODER \
    --ALIGNERS blat,gmap \
    --stringent_alignment_overlap 0.30 \
    --CPU "$cpus"
log "  PASA pipeline complete."

# --- Step 3: Validate PASA output and log assembly counts ---
log "Step 3: Validating PASA output..."
PASA_GFF="${genotype}.pasa_assemblies.gff3"
if [[ -f "$PASA_GFF" ]]; then
    ASSEMBLY_COUNT=$(grep -c '^>' "${PASA_GFF%.gff3}.fasta" 2>/dev/null \
        || grep -c $'\tassembly\t' "$PASA_GFF" || true)
    log "  PASA assemblies: $ASSEMBLY_COUNT"
else
    log "  WARNING: Expected PASA GFF3 not found: $PASA_GFF — check PASA log files."
fi

# --- Summary ---
log "=== PASA pipeline complete ==="
log "  Database:         $DATABASE"
log "  Cleaned FASTA:    $trinity_cleaned ($CLEANED_COUNT transcripts)"
log "  Genome:           $genome"
log "  Aligners used:    BLAT + GMAP"

# =============================================================================
# PASA configuration reference (pasa.config key settings)
# DATABASE=${PWD}/pasa_<genotype>/pasa.sqlite   <- patched into runtime config above
#
# run_spliced_aligners.pl:-N=5
# validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
# validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
# validate_alignments_in_db.dbi:--MAX_INTRON_LENGTH=20000
# subcluster_builder.dbi:-m=50
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=70
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=40
# cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=70
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=50
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=50
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=50
# cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=80
# cDNA_annotation_comparer.dbi:--MAX_UTR_EXON
# =============================================================================
