#!/bin/bash
set -euo pipefail

# Logging helpers
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO  $*" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR $*" >&2; exit 1; }

check_tool() {
    local tool="$1"
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found in PATH: $tool"
    log "Tool check OK: $tool  ($(command -v "$tool"))"
}

# Thread count pulled from SLURM; fall back to 28 for interactive runs
cpus=${SLURM_CPUS_PER_TASK:-28}

# Proteome FASTA files (all must follow the *.protein.fa naming convention)
declare -a SAMPLES=(EV_mazia EV_bedadeti musa_ac musa_ba)
PROT_DIR="EV_Musa_proteins"                              # per-sample protein FASTAs
DB_DIR="db_fasta"                                        # diamond database source FASTAs
OUT_DIR="AnnaRefPlant_NLR.out"                          # main output directory
IPRSC_OUT="interproscan_PF00931_hmmsearch"              # InterProScan output directory
PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_HMM="Pfam-A.hmm"                                   # decompressed HMM flat file
PFAM_NLR_IDS="Pfam_NLR_ID.txt"                         # NLR-related Pfam accessions
PFAM_NB_ARC_ID="PF00931.txt"                            # NB-ARC domain accession
# NLR reference: Kourelis 2021 doi:10.1016/j.molp.2021.08.001; van de Weyer 2019 doi:10.1371/journal.pbio.3001124
NLR_REF_FASTA="${DB_DIR}/ANN_RePlant_NLR.fasta"
INTERPROSCAN_DIR="${INTERPROSCAN_DIR:?INTERPROSCAN_DIR must be set (e.g. /path/to/interproscan-5.xx-xx.0)}"
NLR_ANNOTATOR_MEME="${NLR_ANNOTATOR_MEME:?NLR_ANNOTATOR_MEME must be set (e.g. /path/to/NLR-Annotator/bin/meme.xml)}"

# Module environment
module purge
module load GCC/9.3.0
module load OpenMPI/4.0.3
# Tool availability checks
for tool in wget tar hmmpress hmmfetch hmmbuild hmmsearch diamond clustalo \
            bedtools "${INTERPROSCAN_DIR}/interproscan.sh" fimo mast; do
    check_tool "$tool"
done

# Capture tool versions for reproducibility
log "Tool versions:"
diamond version 2>&1 | head -1 | sed 's/^/  diamond:  /'  >&2 || true
hmmsearch -h 2>&1 | grep -m1 'HMMER'                        >&2 || true
clustalo --version 2>&1 | head -1 | sed 's/^/  clustalo: /' >&2 || true
bedtools --version 2>&1 | head -1 | sed 's/^/  bedtools: /' >&2 || true
log "Using ${cpus} CPUs per SLURM task"
# Input validation
[[ -f "$NLR_REF_FASTA" ]]      || err "NLR reference FASTA not found: $NLR_REF_FASTA"
[[ -f "$PFAM_NLR_IDS" ]]       || err "Pfam NLR ID list not found: $PFAM_NLR_IDS"
[[ -f "$PFAM_NB_ARC_ID" ]]     || err "Pfam NB-ARC ID file not found: $PFAM_NB_ARC_ID"
[[ -f "$NLR_ANNOTATOR_MEME" ]] || err "NLR-Annotator MEME model not found: $NLR_ANNOTATOR_MEME"
for sample in "${SAMPLES[@]}"; do
    fa="${PROT_DIR}/${sample}.protein.fa"
    [[ -f "$fa" && -s "$fa" ]] || err "Protein FASTA missing or empty: $fa"
done

# ── Step 1: Obtain and press Pfam-A HMM database ─────────────────────────────
# NLR Pfam domains: NB-ARC (PF00931), TIR (PF01582, PF13676), RPW8 (PF05659),
# Late blight R1 (PF12061), Rx-type (PF18052), LRR families (PF00560, PF13516,
# PF18805, PF13855, PF07725, PF12799, PF01463, PF01462, PF18831, PF18837,
# PF08263, PF07723, PF13306)
if [[ -f "${PFAM_HMM}.h3m" ]]; then
    log "Pfam-A HMM binary index already present — skipping download and hmmpress"
else
    log "Downloading Pfam-A HMM database from EBI..."
    wget --no-verbose -O "${PFAM_HMM}.gz" "$PFAM_URL"
    gunzip "${PFAM_HMM}.gz"
    hmmpress "$PFAM_HMM"
fi

# Extract NLR-relevant profile HMMs from the full Pfam database
hmmfetch -f "$PFAM_HMM" "$PFAM_NLR_IDS" > Pfam_NLR.hmm
[[ -s Pfam_NLR.hmm ]] || err "hmmfetch produced empty Pfam_NLR.hmm"

# Extract NB-ARC (PF00931) profile separately for second-pass specificity filter
hmmfetch -f "$PFAM_HMM" "$PFAM_NB_ARC_ID" > PF00931.hmm
[[ -s PF00931.hmm ]] || err "hmmfetch produced empty PF00931.hmm"

# ── Step 2: Build diamond database from NLR reference set ────────────────────
log "Building DIAMOND database from annotated NLR reference proteins..."
mkdir -p "$OUT_DIR" "$IPRSC_OUT"
NLR_DB_PREFIX="${NLR_REF_FASTA%.fasta}"   # db prefix without extension → writes .dmnd suffix
diamond makedb --in "$NLR_REF_FASTA" -d "$NLR_DB_PREFIX"
log "DIAMOND database built: ${NLR_DB_PREFIX}.dmnd"

# ── Per-sample NLR prediction loop ───────────────────────────────────────────
declare -A SUMMARY_BLAST SUMMARY_BROAD SUMMARY_NBARC
for sample in "${SAMPLES[@]}"; do
    log "======= Processing sample: ${sample} ======="
    QUERY_FA="${PROT_DIR}/${sample}.protein.fa"
    # 2a. DIAMOND blastp: identify initial NLR candidate proteins
    # --ultra-sensitive: detects divergent NLRs (<30% identity); --masking 0:
    # preserves NLR signal in coiled-coil/LRR regions; -E 1e-5: permissive
    # threshold refined downstream by HMM domain filtering
    log "[${sample}] Running DIAMOND blastp vs. NLR reference..."
    diamond blastp \
        --query   "$QUERY_FA" \
        --db      "$NLR_DB_PREFIX" \
        --ultra-sensitive \
        --masking 0 \
        --out     "${OUT_DIR}/${sample}.blastp.out" \
        --outfmt  6 \
        --compress 0 \
        --evalue  1e-5 \
        --threads "$cpus"
    SUMMARY_BLAST[$sample]=$(wc -l < "${OUT_DIR}/${sample}.blastp.out")
    log "[${sample}] DIAMOND hits: ${SUMMARY_BLAST[$sample]}"

    # 2b. Convert BLAST tabular output to sorted, merged BED coordinates
    # Normalise qstart/qend so start < end; retrieve merged hit subsequences.
    # bedtools getfasta requires a .fai index; create one if not already present.
    awk 'BEGIN{OFS="\t"} {
        s = ($7 < $8) ? $7 : $8
        e = ($7 < $8) ? $8 : $7
        print $1, s-1, e, $3
    }' "${OUT_DIR}/${sample}.blastp.out" \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > "${OUT_DIR}/${sample}.blastp.sorted.merged.bed"
    [[ -f "${QUERY_FA}.fai" ]] || bedtools faidx "$QUERY_FA"
    bedtools getfasta \
        -fi  "$QUERY_FA" \
        -bed "${OUT_DIR}/${sample}.blastp.sorted.merged.bed" \
        > "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.fasta"

    # 2c. Build sample-specific NLR HMM profile from Clustal Omega alignment
    log "[${sample}] Running Clustal Omega alignment for HMM construction..."
    clustalo \
        --in      "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.fasta" \
        --out     "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.aln.fasta" \
        --threads "$cpus" \
        --seqtype Protein
    hmmbuild \
        "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.hmm" \
        "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.aln.fasta"
    # Combine Pfam NLR + sample-specific profile for maximum sensitivity
    cat Pfam_NLR.hmm \
        "${OUT_DIR}/${sample}.blastp.AnnaRefPlant_NLR.hmm" \
        > "${OUT_DIR}/Pfam_RefPlant.${sample}.hmm"

    # 2d. Broad NLR screen with combined HMM library
    # -E 1e-4 / --domE 1e-3: lenient thresholds; domain-envelope coords (cols
    # 20-21) used for boundary extraction — more accurate than alignment coords
    log "[${sample}] Running hmmsearch with combined Pfam+RefPlant NLR profiles..."
    hmmsearch \
        --domtblout "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_dbtblout_e03.out" \
        -E    0.0001 \
        --domE 0.001 \
        --cpu "$cpus" \
        -o    "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam.e03.txt" \
        "${OUT_DIR}/Pfam_RefPlant.${sample}.hmm" \
        "$QUERY_FA"
    grep -v '^#' "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_dbtblout_e03.out" \
        | awk 'BEGIN{OFS="\t"} {print $1, $20-1, $21}' \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed"
    bedtools getfasta \
        -fi  "$QUERY_FA" \
        -bed "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed" \
        > "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_envelope_coordinates.fasta"
    SUMMARY_BROAD[$sample]=$(grep -c '>' "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_envelope_coordinates.fasta" || true)
    log "[${sample}] Broad NLR candidates (post-merge): ${SUMMARY_BROAD[$sample]}"

    # 2e. Second-pass NB-ARC (PF00931) filter — bona fide NLR criterion
    # NB-ARC is universally shared across all NLR subclasses (TIR/CC/RPW8-NLR)
    log "[${sample}] Running NB-ARC (PF00931) specificity filter..."
    hmmsearch \
        --domtblout "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_dbtblout_e03.out" \
        -E    0.0001 \
        --domE 0.001 \
        --cpu "$cpus" \
        -o    "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam.e03.txt" \
        PF00931.hmm \
        "${OUT_DIR}/${sample}.hmmbuild.hmmsearch_pfam_envelope_coordinates.fasta"
    grep -v '^#' "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_dbtblout_e03.out" \
        | awk 'BEGIN{OFS="\t"} {print $1, $20-1, $21}' \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed"
    bedtools getfasta \
        -fi  "$QUERY_FA" \
        -bed "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.sorted.merged.bed" \
        > "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta"
    SUMMARY_NBARC[$sample]=$(grep -c '>' "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta" || true)
    log "[${sample}] High-confidence NLR candidates (NB-ARC confirmed): ${SUMMARY_NBARC[$sample]}"

    # 2f. InterProScan domain annotation of high-confidence NLRs
    # Independent domain evidence (TIR, CC, LRR, RPW8) for NLR subclass classification
    log "[${sample}] Running InterProScan on NB-ARC-confirmed NLRs..."
    "${INTERPROSCAN_DIR}/interproscan.sh" \
        -i   "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta" \
        -f   gff3,TSV \
        -d   "${IPRSC_OUT}/${sample}" \
        --cpu "$cpus"
    # 2g. Motif-based NLR module prediction with NLR-Annotator
    # FIMO: per-motif occurrences; MAST: integrated model scores (doi:10.1094/mpmi-06-18-0165-ta)
    log "[${sample}] Running FIMO + MAST (NLR-Annotator model)..."
    fimo \
        -o "fimo.PF00931.${sample}.out" \
        "$NLR_ANNOTATOR_MEME" \
        "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta"
    mast \
        "$NLR_ANNOTATOR_MEME" \
        "${OUT_DIR}/${sample}.PF00931.hmmsearch_pfam_envelope_coordinates.fasta" \
        -o "mast.${sample}.PF00931.out"

    log "======= Sample ${sample} complete ======="
done

# ── Summary statistics ────────────────────────────────────────────────────────
log "========== NLR PREDICTION SUMMARY =========="
printf "%-15s %12s %12s %12s\n" "Sample" "BLAST_hits" "Broad_NLR" "NB-ARC_NLR" >&2
log "--------------------------------------------"
for sample in "${SAMPLES[@]}"; do
    printf "%-15s %12s %12s %12s\n" \
        "$sample" \
        "${SUMMARY_BLAST[$sample]:-NA}" \
        "${SUMMARY_BROAD[$sample]:-NA}" \
        "${SUMMARY_NBARC[$sample]:-NA}" >&2
done
log "Output directory : ${OUT_DIR}"
log "InterProScan out : ${IPRSC_OUT}"
log "Pipeline complete."
