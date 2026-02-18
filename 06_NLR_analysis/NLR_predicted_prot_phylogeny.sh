#!/bin/bash
# NLR_predicted_prot_phylogeny.sh — NB-ARC domain maximum-likelihood phylogeny.
# Steps: Clustal Omega MSA -> trimAl gap-trim (90% occupancy) -> RAxML (PROTGAMMAAUTO, 1000 bootstraps).
# Input:  EV_MAB.NB_ARC.RPW8_non_plantNBARC.fasta
# Output: RAxML_* files (${RAXML_LABEL}.*)
# Dependencies: clustalo, trimal >= 1.4.1, raxmlHPC-PTHREADS-AVX2
set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO  $*" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR $*" >&2; exit 1; }
check_tool() {
    command -v "$1" >/dev/null 2>&1 || err "Required tool not found in PATH: $1"
    log "Tool OK: $1  ($(command -v "$1"))"
}

cpus=${SLURM_CPUS_PER_TASK:-48}
INPUT_FASTA="EV_MAB.NB_ARC.RPW8_non_plantNBARC.fasta"
ALN_PHY="EV_MAB.NB_ARC.RPW8_non_plantNBARC.phy"
TRIMMED_PHY="EV_MAB.NB_ARC.RPW8_non_plantNBARC_gt_90.phy"
RAXML_LABEL="EV_MAB.NB_ARC.RPW8_non_plantNBARC_90_RAxML"
TRIMAL_BIN="${HOME}/apps/trimal-1.4.1/source/trimal"
RAXML_BIN="raxmlHPC-PTHREADS-AVX2"
RAXML_SEED_X=123456   # bootstrap random seed (-x); fixed for reproducibility
RAXML_SEED_P=246810   # parsimony start-tree seed (-p); fixed for reproducibility
RAXML_BOOTSTRAPS=1000

module purge
module load GCC/9.3.0
module load OpenMPI/4.0.3

check_tool clustalo
[[ -x "$TRIMAL_BIN" ]] || err "trimAl binary not found or not executable: $TRIMAL_BIN"
log "Tool OK: trimal  ($TRIMAL_BIN)"
command -v "$RAXML_BIN" >/dev/null 2>&1 \
    || err "RAxML binary not found: ${RAXML_BIN}. Install PTHREADS-AVX2 build or adjust RAXML_BIN."
log "Tool OK: ${RAXML_BIN}  ($(command -v "$RAXML_BIN"))"

log "Versions:"
clustalo --version 2>&1 | head -1 | sed 's/^/  clustalo: /' >&2 || true
"$TRIMAL_BIN" --version 2>&1 | head -1 | sed 's/^/  trimal:   /' >&2 || true
"$RAXML_BIN" -v 2>&1 | grep -m1 'version' | sed 's/^/  raxml:    /' >&2 || true
log "CPUs: ${cpus}"

[[ -f "$INPUT_FASTA" && -s "$INPUT_FASTA" ]] || err "Input FASTA missing or empty: $INPUT_FASTA"
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA" || true)
log "Input: ${SEQ_COUNT} sequences  (${INPUT_FASTA})"
[[ "$SEQ_COUNT" -ge 4 ]] || err "Fewer than 4 sequences in ${INPUT_FASTA}; RAxML requires >=4 taxa"

# Step 1: Clustal Omega MSA (PHYLIP; 100 HMM-refinement iterations; protein scoring)
log "Running Clustal Omega MSA..."
clustalo \
    --in         "$INPUT_FASTA" \
    --out        "$ALN_PHY" \
    --outfmt=phylip \
    --threads    "$cpus" \
    --iterations 100 \
    --seqtype=Protein
[[ -f "$ALN_PHY" && -s "$ALN_PHY" ]] || err "Clustal Omega produced no output: $ALN_PHY"
log "Alignment: ${ALN_PHY}"

# Step 2: trimAl gap-trim (-gt 0.90 retains columns with >=90% occupancy)
log "Trimming alignment (90% occupancy threshold)..."
"$TRIMAL_BIN" -in "$ALN_PHY" -out "$TRIMMED_PHY" -phylip -gt 0.90
[[ -f "$TRIMMED_PHY" && -s "$TRIMMED_PHY" ]] || err "trimAl produced no output: $TRIMMED_PHY"
TRIMMED_COLS=$(awk 'NR==1{print $2}' "$TRIMMED_PHY" || true)
log "Trimmed: ${TRIMMED_COLS} positions retained  (${TRIMMED_PHY})"

# Step 3: RAxML ML phylogeny (PROTGAMMAAUTO best-fit model; rapid bootstrap -f a)
log "Running RAxML (${RAXML_BOOTSTRAPS} bootstrap replicates)..."
"$RAXML_BIN" \
    -T "$cpus" \
    -s "$TRIMMED_PHY" \
    -n "$RAXML_LABEL" \
    -m PROTGAMMAAUTO \
    -f a \
    -# "$RAXML_BOOTSTRAPS" \
    -x "$RAXML_SEED_X" \
    -p "$RAXML_SEED_P"

BEST_TREE="RAxML_bestTree.${RAXML_LABEL}"
BOOT_FILE="RAxML_bootstrap.${RAXML_LABEL}"
BIPART_TREE="RAxML_bipartitions.${RAXML_LABEL}"
[[ -f "$BEST_TREE" ]]   || err "RAxML produced no best tree: $BEST_TREE"
[[ -f "$BOOT_FILE" ]]   || err "RAxML produced no bootstrap file: $BOOT_FILE"
[[ -f "$BIPART_TREE" ]] || err "RAxML produced no bipartitions tree: $BIPART_TREE"
BOOT_COUNT=$(grep -c '^(' "$BOOT_FILE" || true)

log "=== NLR PHYLOGENY COMPLETE ==="
log "Input          : ${SEQ_COUNT} seqs  (${INPUT_FASTA})"
log "Alignment      : ${ALN_PHY}"
log "Trimmed        : ${TRIMMED_PHY}  (${TRIMMED_COLS:-?} positions)"
log "Best tree      : ${BEST_TREE}"
log "Bootstrap trees: ${BOOT_FILE}  (${BOOT_COUNT} replicates)"
log "Bipartitions   : ${BIPART_TREE}  (use for visualisation)"
log "Model info     : RAxML_info.${RAXML_LABEL}"
