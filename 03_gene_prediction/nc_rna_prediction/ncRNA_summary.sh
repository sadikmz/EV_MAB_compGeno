#!/bin/bash

# ncRNA summary: cascaded non-redundancy pipeline (3 tiers).
# TIER 1 — tRNAscan-SE: gold standard tRNAs, accepted unconditionally.
# TIER 2 — Infernal/cmscan: fills tRNA gaps + all other Rfam ncRNAs (rRNA, snoRNA, snRNA, miRNA…).
# TIER 3 — BLASTN rRNA homology: fills residual rRNA gaps not covered by tiers 1-2.
# Final output per genome: single BED merging all three tiers.
# Usage: bash ncRNA_summary.sh

set -euo pipefail
trap 'rm -f "${TMPDIR:-/tmp}/$$".*' EXIT

# --- Helper functions ---
log()        { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err()        { log "ERROR: $*"; exit 1; }
check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool '$tool' not found in PATH"
    log "  $tool: $("$tool" --version 2>&1 | head -1 || true)"
}
check_file() {
    local filepath=$1 label=$2
    [[ -f "$filepath" ]] || err "Required input not found (${label}): ${filepath}"
    [[ -s "$filepath" ]] || log "WARNING: Input file is empty (${label}): ${filepath}"
}

# --- Configuration ---
GENOTYPES=("mazia" "bedadeti")
TRNASCAN_DIR="tRNAscan"
INFERNAL_DIR="infernal"
HOMOLOGY_DIR="homology_search"
VERSIONS_LOG="${TRNASCAN_DIR}/ncRNA_summary_versions.log"

# --- Pre-flight checks ---
log "=== ncRNA cascaded non-redundancy pipeline starting ==="
for tool in bedtools awk wc; do check_tool "$tool"; done
mkdir -p "$TRNASCAN_DIR" "$INFERNAL_DIR" "$HOMOLOGY_DIR"
{
    echo "=== ncRNA_summary.sh ==="
    bedtools --version 2>&1 | head -1 || true
    echo "Run date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Genotypes: ${GENOTYPES[*]}"
} >> "$VERSIONS_LOG"

# --- Per-genotype processing ---
for genotype in "${GENOTYPES[@]}"; do

    log "--- Processing genotype: ${genotype} ---"

    # Input paths
    TRNASCAN_BED="${TRNASCAN_DIR}/${genotype}.tRNAscan.bed"
    INFERNAL_TBLOUT="${INFERNAL_DIR}/${genotype}.cmscan.tblout"
    INFERNAL_NONTBLOUT="${INFERNAL_DIR}/${genotype}.cmscan.non_tRNA.tblout"
    HOMOLOGY_MERGED_BED="${HOMOLOGY_DIR}/${genotype}.rRNA.e05.blastn.sorted.merged.bed"

    # Output paths
    INFERNAL_TRNA_TBLOUT="${INFERNAL_DIR}/${genotype}.cmscan.tRNA.tblout"
    INFERNAL_TRNA_BED="${INFERNAL_DIR}/${genotype}.cmscan.tRNA.bed"
    INFERNAL_TRNA_SORTED_BED="${INFERNAL_DIR}/${genotype}.cmscan.tRNA.sorted.bed"
    INFERNAL_NON_OL_TRNA_BED="${INFERNAL_DIR}/${genotype}.cmscan.non_overlapping_tRNA.bed"
    HOMOLOGY_NON_OL_BED="${HOMOLOGY_DIR}/${genotype}.rRNA.e05.blastn.non_overlapping.bed"
    HOMOLOGY_NON_OL_INFERNAL_BED="${HOMOLOGY_DIR}/${genotype}.rRNA.e05.blastn.non_overlapping.tRNAscan_infernal.bed"
    FINAL_BED="${TRNASCAN_DIR}/${genotype}_infernal_tRNAscan_tRNA.bed"

    check_file "$TRNASCAN_BED"        "tRNAscan-SE BED"
    check_file "$INFERNAL_TBLOUT"     "Infernal cmscan tblout"
    check_file "$HOMOLOGY_MERGED_BED" "BLASTN rRNA merged BED"

    # TIER 2a — Extract tRNA hits from Infernal tblout and convert to BED
    log "  Extracting Infernal tRNA hits..."
    grep 'tRNA' "$INFERNAL_TBLOUT" > "$INFERNAL_TRNA_TBLOUT" || true   # grep exits 1 if no matches

    # tblout fmt2: seq_from=$10 seq_to=$11 (1-based inclusive). Normalise minus-strand
    # coords (seq_from > seq_to on minus strand) → 0-based half-open BED.
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>0 {
        if ($10 > $11) { s = $11 - 1; e = $10 }
        else           { s = $10 - 1; e = $11 }
        print $4, s, e
    }' "$INFERNAL_TRNA_TBLOUT" > "$INFERNAL_TRNA_BED"
    bedtools sort -i "$INFERNAL_TRNA_BED" > "$INFERNAL_TRNA_SORTED_BED"

    # Extract non-tRNA Rfam hits for tier-3 subtraction and convert to BED
    grep -v '^#\|tRNA' "$INFERNAL_TBLOUT" > "$INFERNAL_NONTBLOUT" || true
    INFERNAL_NONTBLOUT_BED="${INFERNAL_DIR}/${genotype}.cmscan.non_tRNA.bed"
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>0 {
        if ($10 > $11) { s = $11 - 1; e = $10 }
        else           { s = $10 - 1; e = $11 }
        print $4, s, e
    }' "$INFERNAL_NONTBLOUT" \
        | bedtools sort -i - \
        > "$INFERNAL_NONTBLOUT_BED"

    # TIER 2b — Infernal tRNA loci not covered by tRNAscan-SE (gap-fill)
    log "  Subtracting tRNAscan-SE from Infernal tRNA hits..."
    bedtools intersect \
        -a "$INFERNAL_TRNA_SORTED_BED" \
        -b "$TRNASCAN_BED" \
        -v -sorted \
        > "$INFERNAL_NON_OL_TRNA_BED"
    N_INFERNAL_EXTRA=$(wc -l < "$INFERNAL_NON_OL_TRNA_BED")
    log "  Infernal tRNA loci not in tRNAscan-SE (tier-2 additions): ${N_INFERNAL_EXTRA}"

    # TIER 3a — BLASTN rRNA loci not covered by tRNAscan-SE
    log "  Subtracting tRNAscan-SE from BLASTN rRNA hits..."
    bedtools intersect \
        -a "$HOMOLOGY_MERGED_BED" \
        -b "$TRNASCAN_BED" \
        -v -sorted \
        > "$HOMOLOGY_NON_OL_BED"

    # TIER 3b — BLASTN rRNA loci also not covered by Infernal (final gap-fill)
    log "  Subtracting Infernal non-tRNA from residual BLASTN rRNA hits..."
    bedtools intersect \
        -a "$HOMOLOGY_NON_OL_BED" \
        -b "$INFERNAL_NONTBLOUT_BED" \
        -v -sorted \
        > "$HOMOLOGY_NON_OL_INFERNAL_BED"
    N_HOMOLOGY_EXTRA=$(wc -l < "$HOMOLOGY_NON_OL_INFERNAL_BED")
    log "  BLASTN rRNA loci not in tRNAscan-SE or Infernal (tier-3 additions): ${N_HOMOLOGY_EXTRA}"

    # Merge all three tiers into the final non-redundant ncRNA BED
    cat "$TRNASCAN_BED" "$INFERNAL_NON_OL_TRNA_BED" "$HOMOLOGY_NON_OL_INFERNAL_BED" > "$FINAL_BED"
    log "  Final merged ncRNA BED written: ${FINAL_BED}"

done

# --- Summary table ---
log ""
log "=== Final BED file summary ==="
printf "%-12s  %10s  %s\n" "Genotype" "Intervals" "File" >&2
for genotype in "${GENOTYPES[@]}"; do
    FINAL_BED="${TRNASCAN_DIR}/${genotype}_infernal_tRNAscan_tRNA.bed"
    if [[ -f "$FINAL_BED" ]]; then
        N=$(wc -l < "$FINAL_BED")
        printf "%-12s  %10d  %s\n" "$genotype" "$N" "$FINAL_BED" >&2
    else
        printf "%-12s  %10s  %s\n" "$genotype" "MISSING" "$FINAL_BED" >&2
    fi
done

log "=== ncRNA cascaded non-redundancy pipeline complete ==="
