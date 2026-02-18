#!/bin/bash
set -euo pipefail

# NLR-Annotator pipeline: ChopSequence -> NLR-Parser -> NLR-Annotator
# (Steuernagel et al. 2020, Plant Physiology). Dependencies: java>=11, bedtools>=2.29, python3.

# --- Configuration — override via environment variables or edit paths here ---
NLR_PROG="${NLR_PROG:-/path/to/NLR-Annotator}"   # must contain *.jar, meme.xml, meme_4.9.0/
GENOME_DIR="${GENOME_DIR:-/path/to/ref-genomes}"  # expected: mazia.fna bedadeti.fna musa_ac.fna musa_ba.fna
JAVA_MEM="${JAVA_MEM:-32G}"                        # heap for ChopSequence and NLR-Parser (override with JAVA_MEM=Xg)
CPUS="${SLURM_CPUS_PER_TASK:-48}"                 # NLR-Parser MAST worker threads
OUTDIR="${OUTDIR:-NLR-annotator.out}"
IFS=',' read -r -a GENOMES <<< "${GENOMES:-mazia,bedadeti,musa_ac,musa_ba}"
ANNOTATE_SCRIPT="$(dirname "$(realpath "$0")")/annotateClasses.py"
MAST_BIN="${NLR_PROG}/meme_4.9.0/src/mast"

# --- Helper functions ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO  $*" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR $*" >&2; exit 1; }

check_tool() {
    local tool="$1"
    if ! command -v "$tool" &>/dev/null; then err "Required tool not found on PATH: ${tool}"; fi
    log "Tool OK: ${tool} ($(${tool} --version 2>&1 | head -1 || true))"
}

check_jar() {
    local jar="$1"
    if [[ ! -f "$jar" ]]; then err "Required JAR not found: ${jar}"; fi
    log "JAR OK: ${jar}"
}

# --- Input validation ---
log "=== NLR-Annotator pipeline start === CPUS=${CPUS}  JAVA_MEM=${JAVA_MEM}  OUTDIR=${OUTDIR}"

for tool in java bedtools python3; do check_tool "$tool"; done
for jar in ChopSequence NLR-Parser NLR-Annotator; do check_jar "${NLR_PROG}/${jar}.jar"; done

if [[ ! -f "${NLR_PROG}/meme.xml" ]]; then err "MEME motif file not found: ${NLR_PROG}/meme.xml"; fi
if [[ ! -x "$MAST_BIN" ]]; then err "MAST binary not found or not executable: ${MAST_BIN}"; fi
if [[ ! -f "$ANNOTATE_SCRIPT" ]]; then err "annotateClasses.py not found at: ${ANNOTATE_SCRIPT}"; fi

# Validate all genome FASTAs up front (fail fast before any processing).
for genome in "${GENOMES[@]}"; do
    fna="${GENOME_DIR}/${genome}.fna"
    if [[ ! -f "$fna" ]]; then err "Genome FASTA not found: ${fna}"; fi
    log "Genome OK: ${fna}"
done

mkdir -p "${OUTDIR}"
declare -A NLR_COUNTS

# === Per-genome processing loop ===
for genome in "${GENOMES[@]}"; do
    FNA="${GENOME_DIR}/${genome}.fna"
    CHOPSEQ_OUT="${OUTDIR}/${genome}_chopseq_out.fasta"
    PARSER_XML="${OUTDIR}/${genome}.00.nlr.xml"
    ANNOT_TXT="${OUTDIR}/${genome}.nlr.txt"
    ANNOT_GFF="${OUTDIR}/${genome}.NLR-annotator.out.gff"
    ANNOT_BED="${OUTDIR}/${genome}.NLR-annotator.out.bed"
    ANNOT_MOTIFS="${OUTDIR}/${genome}.NLR-annotator.out.motifs.bed"
    ANNOT_FASTA_ALN="${OUTDIR}/${genome}.NLR-annotator.out.nbarkMotifAlignment.fasta"
    NLR_FASTA="${OUTDIR}/${genome}.NLR-annotator.fasta"
    NLR_CLASS="${OUTDIR}/${genome}.NLR_annotator_NLR_class.txt"

    log "--- Processing genome: ${genome} ---"

    # Step 1: ChopSequence — 20 kb overlapping fragments, 5 kb step; skip if already done.
    if [[ -f "$CHOPSEQ_OUT" ]]; then
        log "[${genome}] ChopSequence output exists — skipping (${CHOPSEQ_OUT})"
    else
        log "[${genome}] Step 1: ChopSequence"
        java -Xmx"${JAVA_MEM}" -jar "${NLR_PROG}/ChopSequence.jar" \
            -i "${FNA}" \
            -o "${CHOPSEQ_OUT}"
        log "[${genome}] ChopSequence complete: ${CHOPSEQ_OUT}"
    fi

    # Step 2: NLR-Parser — MAST motif search across fragments (-t worker threads).
    log "[${genome}] Step 2: NLR-Parser"
    java -Xmx"${JAVA_MEM}" -jar "${NLR_PROG}/NLR-Parser.jar" \
        -t "${CPUS}" \
        -y "${MAST_BIN}" \
        -x "${NLR_PROG}/meme.xml" \
        -i "${CHOPSEQ_OUT}" \
        -c "${PARSER_XML}"
    log "[${genome}] NLR-Parser complete: ${PARSER_XML}"

    # Step 3: NLR-Annotator — assemble motif clusters into gene models (GFF3/BED/FASTA outputs).
    log "[${genome}] Step 3: NLR-Annotator"
    java -Xmx"${JAVA_MEM}" -jar "${NLR_PROG}/NLR-Annotator.jar" \
        -i "${PARSER_XML}" \
        -o "${ANNOT_TXT}" \
        -g "${ANNOT_GFF}" \
        -b "${ANNOT_BED}" \
        -m "${ANNOT_MOTIFS}" \
        -a "${ANNOT_FASTA_ALN}"
    log "[${genome}] NLR-Annotator complete"

    # Step 4: Extract NLR locus sequences from genome FASTA using BED coordinates.
    log "[${genome}] Step 4: bedtools getfasta"
    bedtools getfasta \
        -fi "${FNA}" \
        -bed "${ANNOT_BED}" \
        -fo "${NLR_FASTA}"
    log "[${genome}] NLR FASTA written: ${NLR_FASTA}"

    # Step 5: Classify NLR sub-types (CC-NLR, TIR-NLR, RNL, NL, etc.).
    log "[${genome}] Step 5: annotateClasses.py"
    python3 "${ANNOTATE_SCRIPT}" "${ANNOT_TXT}" > "${NLR_CLASS}"
    log "[${genome}] NLR class annotation written: ${NLR_CLASS}"

    NLR_COUNTS["${genome}"]="$( [[ -f "${ANNOT_BED}" ]] && wc -l < "${ANNOT_BED}" || echo 0 )"
    log "[${genome}] NLR loci predicted: ${NLR_COUNTS[${genome}]}"
done

# === Summary table ===
log ""
log "=== NLR-Annotator summary ==="
printf "%-15s  %s\n" "Genome" "NLR_loci" >&2
printf "%-15s  %s\n" "---------------" "--------" >&2
for genome in "${GENOMES[@]}"; do
    printf "%-15s  %s\n" "${genome}" "${NLR_COUNTS[${genome}]}" >&2
done
log "=== Pipeline complete ==="
