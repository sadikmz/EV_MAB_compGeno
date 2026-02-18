#!/bin/bash
# 01_readsCDS_PAV_analysis.sh — WGR → PAV pipeline
# Steps: SRA download → Trim Galore → bwa-mem2 → sort/filter → Picard dedup
#        → samtools flagstat → Qualimap (WG + gene) → bedtools coverage → bigWig
# Edit the Configuration block below, then: bash 01_readsCDS_PAV_analysis.sh

# ── SLURM directives go here (before set -euo pipefail) ──────────────────────

set -euo pipefail

# ── Trap: remove scratch on abnormal exit ─────────────────────────────────────
TMPDIR_SCRATCH=""
cleanup() {
    local exit_code=$?
    [[ -n "$TMPDIR_SCRATCH" && -d "$TMPDIR_SCRATCH" ]] && rm -rf "$TMPDIR_SCRATCH"
    [[ $exit_code -ne 0 ]] && log "Pipeline exited with code ${exit_code}. Check log above."
}
trap cleanup EXIT

# ── Logging helpers ───────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# ── Configuration ─────────────────────────────────────────────────────────────
WORKDIR=$(pwd)
ref_genome="MA"                                                  # reference prefix (MA.fna)
query_genome="EV_mazia"                                          # query prefix (EV_mazia.fna)
PATH_GENOMES="${HOME}/data/genome_dir"                           # genome files directory
PATH_WGS="${WORKDIR}/${query_genome}_trim_out"                   # trimmed-reads output dir
BIOPROJECT="PRJNA344540"                                         # NCBI BioProject for SRA
SRA_LIST_FILE="genotypes_SRA_ID.txt"                            # one SRA accession per line
PICARD_JAR="${HOME}/miniconda3/envs/picard/share/picard-2.27.5-0/picard.jar"
APPSDIR="${HOME}/apps"
BWA_MEM2="${APPSDIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
PICARD_JAVA_MEM="100G"                                           # ~75% of available RAM
QUALIMAP_JAVA_MEM="200G"
TRIM_QUALITY=30                                                  # Phred quality threshold
cpus=${SLURM_CPUS_PER_TASK:-28}

# ── Derived paths ─────────────────────────────────────────────────────────────
ref="${PATH_GENOMES}/${ref_genome}.fna"
gene_bed="${PATH_GENOMES}/${ref_genome}.gene.bed"
output_dir="${WORKDIR}/panEV_${ref_genome}_out"

# ── Tool checks ───────────────────────────────────────────────────────────────
log "=== PAV analysis pipeline: ${query_genome} reads → ${ref_genome} reference ==="
log "Checking required tools..."
check_tool fasterq-dump
check_tool trim_galore
check_tool cutadapt
check_tool samtools
[[ -x "$BWA_MEM2" ]] || err "bwa-mem2 binary not found or not executable: $BWA_MEM2"
log "  bwa-mem2: $("$BWA_MEM2" version 2>&1 | head -1)"
[[ -f "$PICARD_JAR" ]] || err "Picard JAR not found: $PICARD_JAR"
PICARD_VERSION=$(java -jar "$PICARD_JAR" MarkDuplicates --version 2>&1 | head -1 || true)
log "  Picard: $PICARD_VERSION ($PICARD_JAR)"
check_tool qualimap
check_tool bedtools
check_tool bamCoverage

# ── Input validation ──────────────────────────────────────────────────────────
log "Validating inputs..."
[[ -f "$SRA_LIST_FILE" && -s "$SRA_LIST_FILE" ]] || err "SRA list not found or empty: $SRA_LIST_FILE"
[[ -f "$ref"      && -s "$ref"      ]] || err "Reference FASTA not found or empty: $ref"
[[ -f "${ref}.0123" ]] || err "bwa-mem2 index missing for ${ref}. Run 00_prep_data.sh first."
[[ -f "$gene_bed" && -s "$gene_bed" ]] || err "Gene BED not found or empty: $gene_bed. Run 00_prep_data.sh first."

mapfile -t SRA_ARRAY < "$SRA_LIST_FILE"
TOTAL_SAMPLES=${#SRA_ARRAY[@]}
[[ "$TOTAL_SAMPLES" -gt 0 ]] || err "No SRA accessions found in $SRA_LIST_FILE"

log "  Reference genome: $ref"
log "  Gene BED:         $gene_bed ($(wc -l < "$gene_bed") intervals)"
log "  SRA accessions:   $TOTAL_SAMPLES (from $SRA_LIST_FILE)"
log "  Output directory: $output_dir"
log "  CPUs:             $cpus"

mkdir -p "$output_dir" "$PATH_WGS"

PROCESSED=0
SKIPPED_DOWNLOAD=0
FAILED_SAMPLES=()

# ── Per-sample loop ───────────────────────────────────────────────────────────
for READS in "${SRA_ARRAY[@]}"; do

    log "------------------------------------------------------------"
    log "Processing: $READS (sample $((PROCESSED + 1)) of $TOTAL_SAMPLES)"

    Fread="${PATH_WGS}/${READS}_1.fastq.gz"
    Rread="${PATH_WGS}/${READS}_2.fastq.gz"
    out="${ref_genome}_${READS}"
    final_bam="${output_dir}/${out}.allMapped.sorted.markdup.bam"
    flagstat_txt="${output_dir}/${out}.allMapped.sorted.markdup.flagstat.txt"

    # Step 1 — Download + trim (idempotent: skip if trimmed reads already exist)
    if [[ -f "$Fread" && -s "$Fread" && -f "$Rread" && -s "$Rread" ]]; then
        log "  Step 1: Trimmed reads present — skipping download + trim."
        SKIPPED_DOWNLOAD=$((SKIPPED_DOWNLOAD + 1))
    else
        log "  Step 1a: Downloading $READS from BioProject ${BIOPROJECT}..."
        fasterq-dump --split-files --skip-technical --threads "$cpus" --outdir "$PATH_WGS" "$READS"

        # Compress raw FASTQs (prefer pigz, fall back to gzip)
        pigz -p "$cpus" "${PATH_WGS}/${READS}_1.fastq" "${PATH_WGS}/${READS}_2.fastq" \
            2>/dev/null || gzip "${PATH_WGS}/${READS}_1.fastq" "${PATH_WGS}/${READS}_2.fastq"

        [[ -f "${PATH_WGS}/${READS}_1.fastq.gz" && -s "${PATH_WGS}/${READS}_1.fastq.gz" ]] \
            || err "fasterq-dump did not produce R1 for $READS"

        log "  Step 1b: Trimming adapters (Trim Galore, Q${TRIM_QUALITY})..."
        trim_galore \
            --cores 8 --quality "$TRIM_QUALITY" --output_dir "$PATH_WGS" --paired \
            "${PATH_WGS}/${READS}_1.fastq.gz" "${PATH_WGS}/${READS}_2.fastq.gz"

        # Rename Trim Galore output to expected naming convention
        mv "${PATH_WGS}/${READS}_1_val_1.fq.gz" "$Fread"
        mv "${PATH_WGS}/${READS}_2_val_2.fq.gz" "$Rread"
        rm -f "${PATH_WGS}/${READS}_1.fastq.gz" "${PATH_WGS}/${READS}_2.fastq.gz"
    fi

    # Steps 2–4 — Align → sort → dedup → index (skipped if final BAM exists)
    if [[ -f "$final_bam" && -s "$final_bam" ]]; then
        log "  Steps 2–4: Final BAM already present — skipping."
    else
        log "  Step 2: Aligning $READS to ${ref_genome} with bwa-mem2..."
        # bwa-mem2 -M: mark split hits as secondary (Picard compatibility)
        # samtools view -bSF 4: BAM output, discard unmapped reads
        "$BWA_MEM2" mem -t "$cpus" -M "$ref" "$Fread" "$Rread" \
        | samtools view -bSF 4 -@ "$cpus" \
        | samtools sort -@ "$cpus" \
            -T "${output_dir}/tmp_sort_${READS}" \
            -o "${output_dir}/${out}.allMapped.sorted.bam"

        [[ -s "${output_dir}/${out}.allMapped.sorted.bam" ]] \
            || err "bwa-mem2 + sort produced an empty BAM for $READS"

        # Step 3 — Picard MarkDuplicates (REMOVE_DUPLICATES=True for PAV depth accuracy)
        log "  Step 3: Removing duplicates with Picard MarkDuplicates..."
        java -Xmx${PICARD_JAVA_MEM} -jar "$PICARD_JAR" MarkDuplicates \
            INPUT="${output_dir}/${out}.allMapped.sorted.bam" \
            O="$final_bam" \
            M="${output_dir}/${out}.allMapped.sorted.markdup.bam.metrics.txt" \
            REMOVE_DUPLICATES=True \
            VALIDATION_STRINGENCY=LENIENT

        rm -f "${output_dir}/${out}.allMapped.sorted.bam"    # free storage
        [[ -s "$final_bam" ]] || err "Picard MarkDuplicates produced an empty BAM for $READS"

        # Step 4 — Index
        samtools index -@ "$cpus" "$final_bam"
    fi

    # Step 5 — samtools flagstat
    log "  Step 5: Computing mapping statistics (samtools flagstat)..."
    [[ ! -f "$flagstat_txt" ]] && samtools flagstat -@ "$cpus" "$final_bam" > "$flagstat_txt"
    MAPPED_PCT=$(grep "mapped (" "$flagstat_txt" | head -1 | grep -oP '[0-9]+\.[0-9]+(?=%)' || echo "NA")
    log "  Mapped reads: ${MAPPED_PCT}% (see $flagstat_txt)"

    # Step 6 — Qualimap whole-genome QC
    log "  Step 6: Qualimap whole-genome QC..."
    QUALIMAP_WG="${output_dir}/qualimap_${READS}"
    if [[ ! -d "$QUALIMAP_WG" ]]; then
        qualimap bamqc -bam "$final_bam" -outdir "$QUALIMAP_WG" \
            -outfile "${READS}_${ref_genome}.qualimap" \
            -sd -c -nt "$cpus" -outformat PDF:HTML -ip \
            --java-mem-size="$QUALIMAP_JAVA_MEM"
        mv "${QUALIMAP_WG}/genome_results.txt" "${QUALIMAP_WG}/${READS}_genome_results.txt"
        log "  Step 6 complete: ${QUALIMAP_WG}/${READS}_genome_results.txt"
    else
        log "  Step 6: Qualimap whole-genome report already exists — skipping."
    fi

    # Step 7 — Qualimap per-gene-region QC (restricts analysis to gene BED intervals)
    log "  Step 7: Qualimap per-gene-region QC..."
    QUALIMAP_GENE="${output_dir}/qualimap_by_region_${READS}"
    if [[ ! -d "$QUALIMAP_GENE" ]]; then
        qualimap bamqc -bam "$final_bam" -outdir "$QUALIMAP_GENE" \
            -outfile "${READS}_${ref_genome}.qualimap_genes" \
            -sd -c -nt "$cpus" -gff "$gene_bed" \
            -oc "${output_dir}/${READS}_${ref_genome}.qualimap_genes_cov.txt" \
            -os "${output_dir}/${READS}_${ref_genome}.qualimap_gene_summary.txt" \
            -outformat PDF:HTML -ip \
            --java-mem-size="$QUALIMAP_JAVA_MEM"
        mv "${QUALIMAP_GENE}/genome_results.txt" "${QUALIMAP_GENE}/${READS}_genome_results.txt"
        log "  Step 7 complete: ${QUALIMAP_GENE}/${READS}_genome_results.txt"
    else
        log "  Step 7: Qualimap per-gene report already exists — skipping."
    fi

    # Step 8 — bedtools coverage -mean per gene (PAV core metric)
    # Genes with mean depth == 0 across all clade accessions → scored ABSENT
    log "  Step 8: Computing mean per-gene coverage for PAV analysis..."
    GENE_COV="${output_dir}/${READS}.gene_coverage_mean.bed"
    if [[ ! -f "$GENE_COV" || ! -s "$GENE_COV" ]]; then
        bedtools coverage -a "$gene_bed" -b "$final_bam" -mean > "$GENE_COV"
        [[ -s "$GENE_COV" ]] || err "bedtools coverage produced an empty file for $READS: $GENE_COV"
        ZERO_COV=$(awk '$4 == 0 { count++ } END { print count+0 }' "$GENE_COV")
        TOTAL_GENES=$(wc -l < "$GENE_COV")
        log "  Step 8 complete: $GENE_COV | genes: $TOTAL_GENES | zero-cov (absent candidates): $ZERO_COV"
    else
        log "  Step 8: Per-gene coverage file already exists — skipping."
    fi

    # Step 9 — deepTools bamCoverage → RPKM-normalised bigWig for genome browser
    log "  Step 9: Generating bigWig coverage track (deepTools bamCoverage)..."
    BW="${output_dir}/${out}.allMapped.coverage.bw"
    if [[ ! -f "$BW" ]]; then
        bamCoverage -b "$final_bam" --numberOfProcessors "$cpus" --normalizeUsing RPKM -o "$BW"
        [[ -s "$BW" ]] || err "bamCoverage produced an empty bigWig for $READS: $BW"
        log "  Step 9 complete: $BW"
    else
        log "  Step 9: bigWig already exists — skipping."
    fi

    PROCESSED=$((PROCESSED + 1))
    log "  Sample $READS complete."

done

# ── Summary ───────────────────────────────────────────────────────────────────
log "============================================================"
log "=== PAV analysis pipeline complete ==="
log "  Reference genome:      $ref_genome"
log "  Query genome prefix:   $query_genome"
log "  BioProject:            $BIOPROJECT"
log "  Total samples in list: $TOTAL_SAMPLES"
log "  Samples processed:     $PROCESSED"
log "  Downloads skipped:     $SKIPPED_DOWNLOAD (trimmed reads already existed)"
log "  Output directory:      $output_dir"
log ""
log "  Key outputs per sample (READS = SRA accession):"
log "    BAM:            ${output_dir}/${ref_genome}_<READS>.allMapped.sorted.markdup.bam"
log "    Flagstat:       ${output_dir}/${ref_genome}_<READS>.allMapped.sorted.markdup.flagstat.txt"
log "    Dup metrics:    ${output_dir}/${ref_genome}_<READS>.allMapped.sorted.markdup.bam.metrics.txt"
log "    Gene coverage:  ${output_dir}/<READS>.gene_coverage_mean.bed  (PAV input)"
log "    BigWig:         ${output_dir}/${ref_genome}_<READS>.allMapped.coverage.bw"
log "    Qualimap (WG):  ${output_dir}/qualimap_<READS>/"
log "    Qualimap (gene):${output_dir}/qualimap_by_region_<READS>/"
log ""
log "  Next step: run 02_identify_uniqueGenes.sh (or pav_analysis.py) using"
log "  the *.gene_coverage_mean.bed files to build the PAV matrix."
