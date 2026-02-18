#!/bin/bash
# 00_prep_data.sh — Reference genome data preparation for comparative genomics
# Pipeline: faidx -> sizes file -> gene BED -> intergenic BED -> combined BED -> bwa-mem2 index
# Outputs written to PATH_GENOMES; tools required: samtools, gff2bed (bedops), bedtools, bwa-mem2

set -euo pipefail

# --- Logging helpers ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# --- Configuration ---
ref_genome="MA"                              # genome prefix (MA.fna / MA.gff3)
PATH_GENOMES="${HOME}/data/genome_dir"       # directory containing FASTA and GFF3
APPSDIR="${HOME}/apps"
BWA_MEM2="${APPSDIR}/bwa-mem2/bwa-mem2"
cpus=${SLURM_CPUS_PER_TASK:-28}

# --- Derived paths ---
ref="${PATH_GENOMES}/${ref_genome}.fna"
gff="${PATH_GENOMES}/${ref_genome}.gff3"
sizes="${PATH_GENOMES}/${ref_genome}.sizes"
gene_bed="${PATH_GENOMES}/${ref_genome}.gene.bed"
intergenic_bed="${PATH_GENOMES}/${ref_genome}.intergenic.bed"
combined_bed="${PATH_GENOMES}/${ref_genome}.gene_and_intergenic.bed"

# --- Tool checks ---
log "=== Reference genome preparation pipeline ==="
for tool in samtools gff2bed bedtools; do check_tool "$tool"; done
[[ -x "$BWA_MEM2" ]] || err "bwa-mem2 binary not found or not executable: $BWA_MEM2"
log "  bwa-mem2: $("$BWA_MEM2" version 2>&1 | head -1)"

# --- Input validation ---
[[ -f "$ref" && -s "$ref" ]] || err "Reference genome FASTA not found or empty: $ref"
[[ -f "$gff" && -s "$gff" ]] || err "Reference GFF3 not found or empty: $gff"
[[ -w "$PATH_GENOMES" ]]     || err "Output directory is not writable: $PATH_GENOMES"
log "All inputs validated."

# --- Step 1: Index genome with samtools faidx ---
# .fai required for per-chromosome lengths and random-access FASTA in downstream steps.
if [[ -f "${ref}.fai" && "${ref}.fai" -nt "$ref" ]]; then
    log "  .fai index already up-to-date — skipping."
else
    samtools faidx "$ref"
    [[ -s "${ref}.fai" ]] || err "samtools faidx produced an empty index: ${ref}.fai"
    log "  Index written: ${ref}.fai"
fi
SEQ_COUNT=$(wc -l < "${ref}.fai")
log "  Sequences in genome: $SEQ_COUNT"

# --- Step 2: Generate genome sizes file ---
# Two-column TSV (chrom<TAB>length) from .fai cols 1-2; required by bedtools complement.
if [[ -f "$sizes" && "$sizes" -nt "${ref}.fai" ]]; then
    log "  Sizes file already up-to-date — skipping."
else
    cut -f1,2 "${ref}.fai" > "$sizes"
    [[ -s "$sizes" ]] || err "Genome sizes file is empty: $sizes"
    log "  Sizes file written: $sizes"
fi
TOTAL_LENGTH=$(awk '{sum+=$2} END {print sum}' "$sizes")
log "  Total genome length: ${TOTAL_LENGTH} bp across ${SEQ_COUNT} sequences"

# --- Step 3: Extract gene intervals from GFF3 ---
# Select only 'gene' features (col 8) to exclude sub-features; sort for downstream merging.
if [[ -f "$gene_bed" && "$gene_bed" -nt "$gff" ]]; then
    log "  Gene BED already up-to-date — skipping."
else
    gff2bed < "$gff" \
        | awk '$8 == "gene" { print $1, $2, $3 }' OFS='\t' \
        | bedtools sort -i - \
        > "$gene_bed"
    [[ -s "$gene_bed" ]] \
        || err "No gene features extracted from GFF3 (check column 3 feature type): $gff"
    log "  Gene BED written: $gene_bed"
fi
GENE_COUNT=$(wc -l < "$gene_bed")
log "  Gene intervals extracted: $GENE_COUNT"

# --- Step 4: Derive intergenic intervals (complement of gene BED) ---
# bedtools complement returns every interval not covered by genes, bounded by sizes file.
if [[ -f "$intergenic_bed" && "$intergenic_bed" -nt "$gene_bed" ]]; then
    log "  Intergenic BED already up-to-date — skipping."
else
    bedtools complement -i "$gene_bed" -g "$sizes" > "$intergenic_bed"
    [[ -s "$intergenic_bed" ]] \
        || err "bedtools complement produced an empty intergenic BED: $intergenic_bed"
    log "  Intergenic BED written: $intergenic_bed"
fi
INTERGENIC_COUNT=$(wc -l < "$intergenic_bed")
log "  Intergenic intervals: $INTERGENIC_COUNT"

# --- Step 5: Combine and sort gene + intergenic BED ---
# Sorted partition of the full genome used as -a target in PAV analysis.
if [[ -f "$combined_bed" && "$combined_bed" -nt "$intergenic_bed" ]]; then
    log "  Combined BED already up-to-date — skipping."
else
    cat "$gene_bed" "$intergenic_bed" | bedtools sort -i - > "$combined_bed"
    [[ -s "$combined_bed" ]] || err "Combined gene + intergenic BED is empty: $combined_bed"
    log "  Combined BED written: $combined_bed"
fi
COMBINED_COUNT=$(wc -l < "$combined_bed")
COMBINED_BP=$(awk '{sum += ($3-$2)} END {print sum}' "$combined_bed")
log "  Combined BED intervals: $COMBINED_COUNT (${COMBINED_BP} bp total)"
# Sanity-check: combined bp should equal total genome length
if [[ "$COMBINED_BP" -ne "$TOTAL_LENGTH" ]]; then
    log "  WARNING: Combined BED bp (${COMBINED_BP}) != genome size (${TOTAL_LENGTH} bp)."
    log "           This may indicate overlapping gene annotations — run bedtools merge"
    log "           on ${gene_bed} before rerunning if strict partitioning is required."
fi

# --- Step 6: Index genome with bwa-mem2 ---
# Creates FM-index files (.0123, .amb, .ann, .bwt.2bit.64, .pac); re-use if up-to-date.
if [[ -f "${ref}.0123" && "${ref}.0123" -nt "$ref" ]]; then
    log "  bwa-mem2 index already up-to-date — skipping."
else
    "$BWA_MEM2" index "$ref"
    [[ -f "${ref}.0123" ]] \
        || err "bwa-mem2 index sentinel not found after indexing: ${ref}.0123"
    log "  bwa-mem2 index written: ${ref}.{0123,amb,ann,bwt.2bit.64,pac}"
fi

# --- Summary ---
log "=== Reference genome preparation complete ==="
log "  Reference genome:     $ref"
log "  GFF3 annotation:      $gff"
log "  FASTA index:          ${ref}.fai  (${SEQ_COUNT} sequences)"
log "  Genome sizes:         $sizes      (total: ${TOTAL_LENGTH} bp)"
log "  Gene BED:             $gene_bed   (${GENE_COUNT} intervals)"
log "  Intergenic BED:       $intergenic_bed  (${INTERGENIC_COUNT} intervals)"
log "  Combined BED:         $combined_bed    (${COMBINED_COUNT} intervals)"
log "  bwa-mem2 index:       ${ref}.{0123,amb,ann,bwt.2bit.64,pac}"
