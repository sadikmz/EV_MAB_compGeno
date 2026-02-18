#!/bin/bash
#SBATCH --job-name=NLR_loci_vs_genes
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=NLR_loci_vs_genes.%J.out
#SBATCH --error=NLR_loci_vs_genes.%J.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

# NOTE: set -euo pipefail MUST appear after all #SBATCH directives.
set -euo pipefail

# =============================================================================
# NLR loci vs. NLR genes: tblastn-based NLR locus discovery
#
# Maps NB-ARC domain proteins (hmmer/PF00931) to each assembly with tblastn,
# converts HSPs to BED, sorts/merges into discrete loci, subtracts
# NLR-Annotator predictions to isolate candidate novel/divergent NLRs, then
# validates residuals by tblastx against NCBI nt.
#
# .R extension is historical — this is a bash script.
# Usage: sbatch NLR_loci_vs_NLR_genes.R  |or|  bash NLR_loci_vs_NLR_genes.R
# Dependencies: BLAST+ >= 2.10, bedtools >= 2.29
# List files required: genomes_db_list.txt, NLRs_pred_list.txt
# =============================================================================

# --- Paths (edit before submitting) ------------------------------------------
GENOME_DIR=/path/to/ref-genomes   # FASTA files: EV_mazia.fna, EV_bedadeti.fna, Musa_acuminata.fna, Musa_balbisiana.fna
QUERY_DIR=../pred_prot/PF00931_out # per-taxon NB-ARC FASTA files (from NLR_predicted_prot.sh)
NCBI_NT_DB=~/NCBI_NT/BLASTDB/nt   # NCBI nt BLAST db for tblastx validation
NLR_ANNOT_DIR=NLR-annotator.out   # NLR-Annotator output (from NLR-annotator.sh)

CPUS=${SLURM_CPUS_PER_TASK:-48}   # tblastn threads
CPUS_BLASTX=${TBLASTX_CPUS:-16}   # tblastx threads (fewer — NCBI nt is I/O-bound)

BLAST_OUTDIR=hmmNLR_mapping_out
BLASTX_OUTDIR=unmapped_NLR_blastx

# Genome keys → FASTA paths; GENOME_ORDER gives deterministic iteration order.
declare -A GENOME_FASTA=(
    [mazia]="${GENOME_DIR}/EV_mazia.fna"
    [bedadeti]="${GENOME_DIR}/EV_bedadeti.fna"
    [musa_ac]="${GENOME_DIR}/Musa_acuminata.fna"
    [musa_ba]="${GENOME_DIR}/Musa_balbisiana.fna"
)
GENOME_ORDER=(mazia bedadeti musa_ac musa_ba)

# --- Helper functions --------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO  $*" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR $*" >&2; }

check_tool() {
    local tool="$1"
    command -v "$tool" &>/dev/null \
        || { err "Required tool not found on PATH: ${tool}"; exit 1; }
    log "Tool OK: ${tool} ($(${tool} -version 2>&1 | head -1 || true))"
}

count_bed() { [[ -f "$1" ]] && wc -l < "$1" || echo 0; }

# --- Input validation --------------------------------------------------------
log "=== NLR loci vs. NLR genes pipeline start === CPUS=${CPUS}  CPUS_BLASTX=${CPUS_BLASTX}"

for tool in makeblastdb tblastn tblastx bedtools awk; do check_tool "$tool"; done

for listfile in genomes_db_list.txt NLRs_pred_list.txt; do
    [[ -f "$listfile" ]] || { err "Required list file not found: ${listfile}"; exit 1; }
    log "List file OK: ${listfile}"
done

for genome in "${GENOME_ORDER[@]}"; do
    fna="${GENOME_FASTA[${genome}]}"
    [[ -f "$fna" ]] || { err "Genome FASTA not found: ${fna} (key: ${genome})"; exit 1; }
    annot_bed="${NLR_ANNOT_DIR}/${genome}.NLR-annotator.out.bed"
    [[ -f "$annot_bed" ]] || { err "NLR-Annotator BED not found: ${annot_bed} — run NLR-annotator.sh first"; exit 1; }
done

mkdir -p "${BLAST_OUTDIR}" "${BLASTX_OUTDIR}" NLR_mapping_out_multihits

# =============================================================================
# Part 1: Build BLAST nucleotide databases
# -parse_seqids preserves sequence identifiers for downstream getfasta steps.
# =============================================================================
log "=== Part 1: Building BLAST nucleotide databases ==="

for genome in "${GENOME_ORDER[@]}"; do
    fna="${GENOME_FASTA[${genome}]}"
    log "makeblastdb: ${fna}"
    makeblastdb -in "${fna}" -dbtype nucl -parse_seqids
done

# =============================================================================
# Part 2: tblastn — map NB-ARC proteins to each genome
#
# tblastn queries protein vs. 6-frame translated genome (more sensitive than
# blastn for divergent NLRs). -max_hsps 1 keeps the best HSP per pair.
# awk normalises minus-strand coordinates (col9 > col10) to BED start < end
# and converts from 1-based to 0-based half-open intervals (subtract 1 from
# start). bedtools sort+merge collapses overlapping HSPs into single loci.
# =============================================================================
log "=== Part 2: tblastn NB-ARC protein mapping ==="

while IFS= read -r genome_key; do
    fna="${GENOME_FASTA[${genome_key}]:-}"
    [[ -n "$fna" ]] || { err "Genome key '${genome_key}' not in GENOME_FASTA map"; exit 1; }

    while IFS= read -r query_prefix; do
        log "  tblastn: query=${query_prefix}  genome=${genome_key}"

        QUERY_FASTA="${QUERY_DIR}/${query_prefix}.NB_ARC.fasta"
        [[ -f "$QUERY_FASTA" ]] || { err "NB-ARC query FASTA not found: ${QUERY_FASTA}"; exit 1; }

        RAW_OUT="${BLAST_OUTDIR}/${query_prefix}.${genome_key}.e5.tblastn.out"
        BED_MERGED="${BLAST_OUTDIR}/${query_prefix}.${genome_key}.NLR_genes.tblastn.sorted.merged.bed"

        tblastn \
            -db "${fna}" -query "${QUERY_FASTA}" \
            -evalue 1e-5 -outfmt 6 -max_hsps 1 \
            -num_threads "${CPUS}" -out "${RAW_OUT}"

        # Swap minus-strand coords → 0-based BED, then sort+merge into loci.
        awk 'BEGIN{OFS="\t"}
             $9 > $10 { print $2, $10-1, $9,  $3; next }
                      { print $2, $9-1,  $10, $3 }' "${RAW_OUT}" \
            | bedtools sort -i - \
            | bedtools merge -i - \
            > "${BED_MERGED}"

        log "  Merged NLR loci: query=${query_prefix}  genome=${genome_key}  loci=$(count_bed "${BED_MERGED}")"

    done < NLRs_pred_list.txt
done < genomes_db_list.txt

# =============================================================================
# Part 3: Identify NLR loci NOT overlapping NLR-Annotator gene models
#
# Union of all per-query tblastn BED files → subtract NLR-Annotator loci
# (bedtools intersect -v = set difference A \ B). Residuals are candidate
# divergent NLRs where the NB-ARC domain was detected but the full motif
# pattern was too degenerate for NLR-Annotator.
# =============================================================================
log "=== Part 3: Non-overlapping NLR locus extraction ==="

declare -A NONOLAP_COUNTS

for genome in "${GENOME_ORDER[@]}"; do
    FNA="${GENOME_FASTA[${genome}]}"
    ANNOT_BED="${NLR_ANNOT_DIR}/${genome}.NLR-annotator.out.bed"
    UNION_BED="${BLAST_OUTDIR}/${genome}.all_queries.union.bed"
    NONOLAP_BED="${BLAST_OUTDIR}/${genome}_NLR_loci_non_overlapping.bed"
    NONOLAP_FASTA="${BLAST_OUTDIR}/${genome}.NLR_loci_non_overlapping.fasta"

    # Merge all per-query loci for this genome into a cross-query union.
    # nullglob: expand to empty array (not literal string) if no files match.
    shopt -s nullglob
    bed_files=( "${BLAST_OUTDIR}/"*".${genome}.NLR_genes.tblastn.sorted.merged.bed" )
    shopt -u nullglob
    [[ ${#bed_files[@]} -gt 0 ]] \
        || err "[${genome}] No tblastn BED files found in ${BLAST_OUTDIR}/ — run Part 2 first"
    cat "${bed_files[@]}" \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > "${UNION_BED}"
    log "  [${genome}] Union NLR loci: $(count_bed "${UNION_BED}")  →  subtracting NLR-Annotator"

    # Set difference: tblastn loci with no NLR-Annotator overlap.
    bedtools intersect -a "${UNION_BED}" -b "${ANNOT_BED}" -v > "${NONOLAP_BED}"
    NONOLAP_COUNTS["${genome}"]=$(count_bed "${NONOLAP_BED}")
    log "  [${genome}] Non-overlapping NLR loci: ${NONOLAP_COUNTS[${genome}]}"

    # Extract sequences for tblastx validation.
    bedtools getfasta -fi "${FNA}" -bed "${NONOLAP_BED}" -fo "${NONOLAP_FASTA}"
done

# =============================================================================
# Part 4: tblastx validation against NCBI nt
#
# Translates both query and subject in all 6 frames — sensitive enough to
# detect distant NLR homologs or TE fragments. Hits to NLR proteins in nt
# confirm loci as genuine NLR remnants rather than assembly artefacts.
# =============================================================================
log "=== Part 4: tblastx validation of non-overlapping NLR loci ==="

for genome in "${GENOME_ORDER[@]}"; do
    QUERY_FASTA="${BLAST_OUTDIR}/${genome}.NLR_loci_non_overlapping.fasta"
    BLASTX_OUT="${BLASTX_OUTDIR}/${genome}.unmapped_NLR.tblastx.e5.out"

    if [[ ! -s "$QUERY_FASTA" ]]; then
        log "  [${genome}] No non-overlapping loci — skipping tblastx"
        continue
    fi

    tblastx \
        -db "${NCBI_NT_DB}" -query "${QUERY_FASTA}" \
        -evalue 1e-5 -outfmt 6 -max_hsps 1 \
        -num_threads "${CPUS_BLASTX}" -out "${BLASTX_OUT}"

    log "  [${genome}] tblastx hits: $(wc -l < "${BLASTX_OUT}" 2>/dev/null || echo 0)"
done

# =============================================================================
# Summary
# =============================================================================
log ""
log "=== NLR loci vs. NLR genes — summary ==="
printf "%-15s  %s\n" "Genome" "Non-overlapping_loci" >&2
printf "%-15s  %s\n" "---------------" "--------------------" >&2
for genome in "${GENOME_ORDER[@]}"; do
    printf "%-15s  %s\n" "${genome}" "${NONOLAP_COUNTS[${genome}]:-0}" >&2
done
log "=== Pipeline complete ==="
