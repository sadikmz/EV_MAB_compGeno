#!/bin/bash
# braker.sh — BRAKER2 ETP gene prediction (ProtHint protein hints + RNA-seq BAM).
# Pipeline: (1) ProtHint hints, (2) env/input validation, (3) BRAKER2 ETP,
#           (4) count gene models, (5) note: copy Augustus species model for MAKER.
# Tools: prothint.py, braker.pl, Augustus, GeneMark-ET
# Proteins: 13 monocot/dicot spp. from Phytozome v2 (Arabidopsis, Oryza, Sorghum,
#   Dioscorea, Brachypodium, Manihot, Zea, Ananas, Musa acuminata/balbisiana/
#   schizocarpa, Ensete). Refs: github.com/Gaius-Augustus/BRAKER, gatech-genemark/ProtHint

set -euo pipefail

# Logging helpers
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# Configuration — edit before running
genotype="${GENOTYPE:?GENOTYPE env var must be set (e.g. export GENOTYPE=ensete_ventricosum)}"
softmasked_genome="softmasked.fna"    # Repeat-softmasked reference genome FASTA
PROTHINT_BIN="${HOME}/apps/ProtHint/bin/prothint.py"
PROTEIN_DB="combined_multi_species_proteins.fa"
BASENAME="${softmasked_genome%.fna}"; RNASEQ_BAM="data/${BASENAME}.rnaseq.bam"
PROTHINT_GFF="prothint_augustus.gff"
BRAKER_WORKDIR="braker"; AUGUSTUS_SPECIES="${genotype}"
cpus=${SLURM_CPUS_PER_TASK:-28}

# Tool checks
log "=== BRAKER2 gene prediction pipeline ==="
check_tool braker.pl
[[ -x "$PROTHINT_BIN" ]] || err "ProtHint not found at: $PROTHINT_BIN"
log "  prothint.py: $("$PROTHINT_BIN" --version 2>&1 | head -1)"

# Environment validation
log "Validating BRAKER2 environment..."
# AUGUSTUS_CONFIG_PATH must be set and writable; Augustus writes trained params there.
[[ -n "${AUGUSTUS_CONFIG_PATH:-}" ]] \
    || err "AUGUSTUS_CONFIG_PATH is not set. Export it before running: export AUGUSTUS_CONFIG_PATH=/path/to/Augustus/config"
[[ -d "$AUGUSTUS_CONFIG_PATH" ]] \
    || err "AUGUSTUS_CONFIG_PATH directory does not exist: $AUGUSTUS_CONFIG_PATH"
[[ -w "$AUGUSTUS_CONFIG_PATH/species" ]] \
    || err "AUGUSTUS_CONFIG_PATH/species is not writable: ${AUGUSTUS_CONFIG_PATH}/species. BRAKER cannot write trained parameters."
log "  AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH (writable)"
# GeneMark requires a licence key in the home directory.
[[ -f "${HOME}/.gm_key" ]] \
    || err "GeneMark licence key not found: ~/.gm_key. Download from http://topaz.gatech.edu/GeneMark/license_download.cgi"
log "  GeneMark licence key: ~/.gm_key found."

# Input validation
log "Validating input files..."
[[ -f "$softmasked_genome" && -s "$softmasked_genome" ]] || err "Softmasked genome not found or empty: $softmasked_genome"
[[ -f "$PROTEIN_DB"        && -s "$PROTEIN_DB"        ]] || err "Protein database not found or empty: $PROTEIN_DB"
[[ -f "$RNASEQ_BAM"        && -s "$RNASEQ_BAM"        ]] || err "RNA-seq BAM not found or empty: $RNASEQ_BAM"
log "All input files validated."

# Step 1 — Generate protein hints with ProtHint
# Aligns protein DB to genome (DIAMOND + Spaln) → Augustus-compatible hint GFF.
log "Step 1: Generating protein hints with ProtHint..."
if [[ -f "$PROTHINT_GFF" && -s "$PROTHINT_GFF" ]]; then
    log "  ProtHint GFF already exists — skipping: $PROTHINT_GFF"
else
    "$PROTHINT_BIN" "$softmasked_genome" "$PROTEIN_DB" --threads "$cpus"
    [[ -f "$PROTHINT_GFF" && -s "$PROTHINT_GFF" ]] \
        || err "ProtHint did not produce a non-empty GFF: $PROTHINT_GFF"
    log "  ProtHint GFF produced: $PROTHINT_GFF"
fi
HINT_COUNT=$(grep -c $'\tCDS\t' "$PROTHINT_GFF" || true)
log "  CDS hints in ProtHint GFF: $HINT_COUNT"

# Step 2 — Run BRAKER2 in ETP mode
# --bam RNA-seq BAM, --hints protein GFF, --etpmode trains Augustus+GeneMark jointly,
# --softmasking informs Augustus of lowercase repeats,
# --AUGUSTUS_ab_initio outputs pure ab initio predictions alongside evidence-based ones.
log "Step 2: Running BRAKER2 (ETP mode)..."
mkdir -p "$BRAKER_WORKDIR"
braker.pl \
    --species="$AUGUSTUS_SPECIES" \
    --genome="$softmasked_genome" \
    --hints="$PROTHINT_GFF" \
    --bam="$RNASEQ_BAM" \
    --etpmode \
    --softmasking \
    --AUGUSTUS_ab_initio \
    --workingdir="$BRAKER_WORKDIR" \
    --threads "$cpus"
log "  BRAKER2 run complete."

# Step 3 — Validate output and log gene model counts
log "Step 3: Validating BRAKER2 output..."
BRAKER_GFF="${BRAKER_WORKDIR}/augustus.hints.gff3"
[[ -f "$BRAKER_GFF" ]] || err "BRAKER2 output GFF3 not found: $BRAKER_GFF"
[[ -s "$BRAKER_GFF" ]] || err "BRAKER2 output GFF3 is empty: $BRAKER_GFF"
GENE_COUNT=$(grep -c $'\tgene\t' "$BRAKER_GFF" || true)
MRNA_COUNT=$(grep -c $'\tmRNA\t' "$BRAKER_GFF" || true)
log "  Gene models in ${BRAKER_GFF}:"
log "    gene features: $GENE_COUNT"
log "    mRNA features: $MRNA_COUNT"
# Also report ab initio prediction count if produced
ABINITIO_GFF="${BRAKER_WORKDIR}/augustus.ab_initio.gff3"
if [[ -f "$ABINITIO_GFF" ]]; then
    ABINITIO_GENE_COUNT=$(grep -c $'\tgene\t' "$ABINITIO_GFF" || true)
    log "  Ab initio gene models: $ABINITIO_GENE_COUNT (${ABINITIO_GFF})"
fi

# Step 4 — Copy trained Augustus species model to shared config
# Uncomment once BRAKER2 finishes. MAKER needs species params at
# AUGUSTUS_CONFIG_PATH/species/<species_name>/ to call --species=<species_name>.
#   TRAINED_SPECIES_DIR="${BRAKER_WORKDIR}/augustus_output/species/${AUGUSTUS_SPECIES}"
#   DEST="${AUGUSTUS_CONFIG_PATH}/species/${AUGUSTUS_SPECIES}"
#   if [[ -d "$TRAINED_SPECIES_DIR" && ! -d "$DEST" ]]; then
#       cp -r "$TRAINED_SPECIES_DIR" "$DEST"
#       log "Copied trained species model to: $DEST"
#   else
#       log "Species model already present at $DEST — skipping copy."
#   fi
log "NOTE: After confirming BRAKER2 output, copy the trained Augustus species"
log "  model from ${BRAKER_WORKDIR}/augustus_output/species/${AUGUSTUS_SPECIES}"
log "  to \${AUGUSTUS_CONFIG_PATH}/species/${AUGUSTUS_SPECIES}"
log "  so that MAKER can locate it via --species=${AUGUSTUS_SPECIES}."

# Summary
log "=== BRAKER2 pipeline complete ==="
log "  Working dir:     $BRAKER_WORKDIR"
log "  Hints GFF:       ${BRAKER_WORKDIR}/augustus.hints.gff3 ($GENE_COUNT genes)"
log "  GeneMark model:  ${BRAKER_WORKDIR}/GeneMark_hmm_full.mod"
log "  ProtHint GFF:    $PROTHINT_GFF"
