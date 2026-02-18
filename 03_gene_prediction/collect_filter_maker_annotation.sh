#!/bin/bash

set -euo pipefail

# Collect MAKER output, merge GFF3/FASTA, filter by AED, assess quality with BUSCO (AED sweep 0.20-0.35).

# --- Parameters (override via environment or edit here) ---
genotype="${GENOTYPE:-genotype_prefix}"
MAKER_OUTPUT_DIR="${MAKER_OUTPUT_DIR:-${genotype}.maker.output}"
MASTER_INDEX="${MASTER_INDEX:-${MAKER_OUTPUT_DIR}/${genotype}_master_datastore_index.log}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-embryophyta_odb10}"
BUSCO_CPU="${BUSCO_CPU:-48}"

# AED sweep: generate one FASTA per threshold from 0.20 to 0.35 inclusive
AED_THRESHOLDS=(0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35)

# --- Logging helpers ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }
check_tool() {
    command -v "$1" >/dev/null 2>&1 || err "Required tool not found: $1"
    log "  $1: $("$1" --version 2>&1 | head -1)"
}

# --- Pre-flight checks ---
log "=== collect_filter_maker_annotation ==="
log "Genotype     : ${genotype}"
log "MAKER index  : ${MASTER_INDEX}"
log "BUSCO lineage: ${BUSCO_LINEAGE}"
for tool in gff3_merge fasta_merge busco generate_plot.py; do check_tool "${tool}"; done
[[ -f "${MASTER_INDEX}" ]] || err "MAKER master datastore index not found: ${MASTER_INDEX}"

# --- Step 1: Merge MAKER outputs into unified GFF3 and FASTA ---
log "--- Merging MAKER outputs ---"
mkdir -p gene_models predicted_fasta gene_utr_count AED_fasta collect_split

gff3_merge \
    -d "${MASTER_INDEX}" \
    -o "gene_models/${genotype}.maker.gff"

fasta_merge \
    -d "${MASTER_INDEX}" \
    -o "predicted_fasta/${genotype}"

# Extract MAKER gene models only (source == "maker"); -g is not a valid gff3_merge flag.
awk '$2 == "maker"' \
    "gene_models/${genotype}.maker.gff" \
    > "gene_models/${genotype}.maker.gene_models.gff"

all_proteins="predicted_fasta/${genotype}.all.maker.proteins.fasta"
[[ -f "${all_proteins}" ]] || err "Merged protein FASTA not found: ${all_proteins}"
total_proteins=$(grep -c '^>' "${all_proteins}" 2>/dev/null || echo 0)
log "  Total predicted proteins: ${total_proteins}"

# --- Step 2: Count genes and UTRs across a coarse AED sweep (GFF-based) ---
log "--- Counting genes/UTRs per AED threshold (GFF-based) ---"
for aed_thresh in 0.16 0.21 0.26 0.31 0.36 0.41 0.46 0.51 0.56 0.61 0.66 0.71 0.76 0.81 0.86 0.91 0.96 1.00; do
    quality_filter.pl -a "${aed_thresh}" "gene_models/${genotype}.maker.gene_models.gff" \
        > "gene_models/maker.${aed_thresh}.gff" 2>/dev/null || true
    gff_file="gene_models/maker.${aed_thresh}.gff"
    if [[ -f "${gff_file}" ]]; then
        grep -cP '\tgene\t'           "${gff_file}" >> gene_utr_count/genes.AED.txt || echo 0 >> gene_utr_count/genes.AED.txt
        grep -cP '\tthree_prime_UTR\t' "${gff_file}" >> gene_utr_count/utr3.AED.txt  || echo 0 >> gene_utr_count/utr3.AED.txt
        grep -cP '\tfive_prime_UTR\t'  "${gff_file}" >> gene_utr_count/utr5.AED.txt  || echo 0 >> gene_utr_count/utr5.AED.txt
    fi
done

paste gene_utr_count/genes.AED.txt gene_utr_count/utr3.AED.txt gene_utr_count/utr5.AED.txt \
    > gene_utr_count/genes.35utr.txt
log "  Gene/UTR count table: gene_utr_count/genes.35utr.txt"

# Generate AED CDF (optional tool — guard to avoid hard failure if not installed)
if command -v AED_cdf_generator.pl &>/dev/null; then
    AED_cdf_generator.pl -b 0.01 \
        "gene_models/${genotype}.maker.gene_models.gff" \
        > AED_unfiltered.genemodels.txt
    log "  AED CDF: AED_unfiltered.genemodels.txt"
else
    log "  WARNING: AED_cdf_generator.pl not found — run manually: AED_cdf_generator.pl -b 0.01 gene_models/${genotype}.maker.gene_models.gff"
fi

# --- Step 3: Filter protein FASTA by AED threshold (sweep 0.20-0.35) ---
# Uses awk for accurate numeric comparison; grep/regex mis-classifies values
# (e.g. character class [0-2] is not a numeric comparison).
log "--- Filtering proteins by AED threshold (awk numeric comparison) ---"

filter_by_aed() {
    local input_fasta="$1" output_fasta="$2" max_aed="$3"
    awk -v threshold="${max_aed}" '
        /^>/ {
            keep = 0
            if (match($0, /AED:([0-9.]+)/, arr)) {
                keep = (arr[1] + 0 <= threshold + 0) ? 1 : 0
            }
            if (keep) print
            next
        }
        keep { print }
    ' "${input_fasta}" > "${output_fasta}"
    local n
    n=$(grep -c '^>' "${output_fasta}" 2>/dev/null || echo 0)
    log "  AED<=${max_aed}: ${n} proteins -> $(basename "${output_fasta}")"
}

for aed in "${AED_THRESHOLDS[@]}"; do
    aed_label=$(printf "le%.2f" "${aed}")
    filter_by_aed \
        "${all_proteins}" \
        "AED_fasta/${genotype}.maker.proteins.AED_${aed_label}.fasta" \
        "${aed}"
done
log "  All AED-filtered FASTAs written to: AED_fasta/"

# --- Step 4: BUSCO quality assessment across AED-filtered FASTAs ---
log "--- Running BUSCO across AED-filtered FASTAs ---"
for aed in "${AED_THRESHOLDS[@]}"; do
    aed_label=$(printf "le%.2f" "${aed}")
    fasta="AED_fasta/${genotype}.maker.proteins.AED_${aed_label}.fasta"
    if [[ ! -f "${fasta}" ]]; then
        log "  WARNING: FASTA not found for AED ${aed_label} -- skipping BUSCO"
        continue
    fi
    n_prot=$(grep -c '^>' "${fasta}" 2>/dev/null || echo 0)
    if [[ "${n_prot}" -eq 0 ]]; then
        log "  WARNING: Empty FASTA for AED ${aed_label} -- skipping BUSCO"
        continue
    fi
    busco_out="busco.AED_${aed_label}.out"
    log "  BUSCO for AED<=${aed}: ${n_prot} proteins -> ${busco_out}"
    busco \
        -i  "${fasta}" \
        -m  protein \
        -l  "${BUSCO_LINEAGE}" \
        -c  "${BUSCO_CPU}" \
        -o  "${busco_out}" \
        --long \
        -f \
        --offline
    summary_src="${busco_out}/short_summary.specific.${BUSCO_LINEAGE}.${busco_out}.txt"
    if [[ -f "${summary_src}" ]]; then
        cp "${summary_src}" "collect_split/"
        log "  Copied BUSCO summary: $(basename "${summary_src}")"
    else
        log "  WARNING: BUSCO summary not found: ${summary_src}"
    fi
done

# --- Step 5: Generate BUSCO comparison plot ---
log "--- Generating BUSCO comparison plot ---"
if [[ -d collect_split ]] && ls collect_split/short_summary.*.txt >/dev/null 2>&1; then
    generate_plot.py -wd collect_split
    log "  BUSCO plot generated in: collect_split/"
else
    log "  WARNING: No BUSCO summaries found in collect_split/ -- skipping plot"
fi

log "=== collect_filter_maker_annotation complete ==="
