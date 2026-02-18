#!/bin/bash
# maker.sh — Iterative genome annotation with MAKER2 (MPI-enabled), 4 rounds:
#   R1: Bootstrap from EST+protein evidence (est2genome=1, protein2genome=1).
#   R2-R4: Retrain SNAP each round; evidence-only mode OFF (est2genome=0, protein2genome=0).
#   R4 (final): SNAP R3 + GeneMark (BRAKER2) + Augustus; merge and functionally annotate.
# Prerequisites: maker_opts.ctl, maker_exe.ctl, maker_bopts.ctl pre-generated with `maker -CTL`.

set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
err() { log "ERROR: $*"; exit 1; }

check_tool() {
    local tool=$1
    command -v "$tool" >/dev/null 2>&1 || err "Required tool not found: $tool"
    log "  $tool: $("$tool" --version 2>&1 | head -1)"
}

# Configuration — override via environment variables or edit here
genotype="${GENOTYPE:-genotype_prefix}"  # Sample / accession identifier
MAKERDIR="${genotype}"                   # MAKER base name; output goes to ${MAKERDIR}.maker.output/
maxIntronLen="${MAX_INTRON_LEN:-82825}"  # Max intron length (bp), empirically determined from RNA-seq
cpus="${SLURM_CPUS_PER_TASK:-28}"
export OMP_NUM_THREADS=1   # MPI parallelism: each rank must not spawn extra threads (cpus² oversubscription)
MAKER_EXE="${MAKER_EXE:-maker_exe.ctl}"
MAKER_OPTS="${MAKER_OPTS:-maker_opts.ctl}"
MAKER_BOPTS="${MAKER_BOPTS:-maker_bopts.ctl}"

# Validate MAKER master datastore index exists and is non-empty
check_maker_output() {
    local round_label="$1"
    local index="${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log"
    [[ -f "$index" ]] || err "[${round_label}] MAKER master datastore index not found: $index"
    [[ -s "$index" ]] || err "[${round_label}] MAKER master datastore index is empty: $index"
    log "  [${round_label}] Datastore index OK — $(wc -l < "$index") scaffold entries."
}

# Tool and control-file checks
for tool in maker mpiexec gff3_merge fasta_merge maker2zff fathom forge hmm-assembler.pl; do
    check_tool "$tool"
done
AED_CDF=$(command -v AED_cdf_generator.pl 2>/dev/null || true)

for ctl in "$MAKER_EXE" "$MAKER_OPTS" "$MAKER_BOPTS"; do
    [[ -f "$ctl" ]] || err "MAKER control file not found: $ctl. Run \`maker -CTL\` first."
done

# ===========================================================================
# ROUND 1 — Bootstrap gene predictions from EST and protein evidence
# maker_opts.ctl: est2genome=1, protein2genome=1, cpu=1
# ===========================================================================
log "=== Round 1: Bootstrap from EST + protein evidence ==="
mkdir -p round1

mpiexec -n "$cpus" maker "$MAKER_EXE" "$MAKER_OPTS" "$MAKER_BOPTS" \
    -qq -base "$MAKERDIR" -fix_nucleotides

check_maker_output "Round 1"

gff3_merge \
    -d "${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log" \
    -o "round1/${genotype}.snap01.all.gff"

[[ -f "round1/${genotype}.snap01.all.gff" ]] \
    || err "Round 1 merged GFF not produced: round1/${genotype}.snap01.all.gff"
log "  Round 1 gene models: $(grep -c $'\tgene\t' "round1/${genotype}.snap01.all.gff" || true)"

# ===========================================================================
# SNAP Training — Round 1
# maker2zff converts MAKER GFF3 to ZFF; fathom validates; forge+hmm-assembler builds HMM
# ===========================================================================
log "=== SNAP training on Round 1 predictions ==="
mkdir -p snap.round1 train_snap.round1

(
    cd snap.round1
    maker2zff "../round1/${genotype}.snap01.all.gff"
    [[ -s genome.ann && -s genome.dna ]] || { echo "ERROR: maker2zff did not produce genome.ann/genome.dna" >&2; exit 1; }
    fathom genome.ann genome.dna -validate > snap_validate_output.txt
    grep "error" snap_validate_output.txt > snap.error.out || true
    grep 'errors' snap_validate_output.txt | awk '{print $1, $2}' > snap.scaffolds_gene_model_error.txt || true
    grep 'errors' snap_validate_output.txt | awk '{print $2}' > snap.gene_model_error.txt || true
    grep -vwE -f snap.gene_model_error.txt genome.ann > genome.ann1 || cp genome.ann genome.ann1
    fathom genome.ann1 genome.dna -validate | grep "error" > snap.summary.txt || true
)

(
    cd train_snap.round1
    fathom ../snap.round1/genome.ann1 ../snap.round1/genome.dna -categorize 1000
    fathom uni.ann uni.dna -export 1000 -plus
    forge export.ann export.dna
    hmm-assembler.pl "${genotype}.snap1" . > "${genotype}.snap1.hmm"
)

[[ -f "train_snap.round1/${genotype}.snap1.hmm" ]] || err "SNAP Round 1 HMM not produced."
log "  SNAP Round 1 HMM: train_snap.round1/${genotype}.snap1.hmm"

# ===========================================================================
# ROUND 2 — SNAP (Round 1 HMM); evidence-only mode OFF
# maker_opts.ctl: est2genome=0, protein2genome=0,
#   maker_gff=<abs>/round1/<genotype>.snap01.all.gff,
#   snaphmm=<abs>/train_snap.round1/<genotype>.snap1.hmm
# ===========================================================================
log "=== Round 2: SNAP Round 1 HMM ==="
mkdir -p round2

mpiexec -n "$cpus" maker "$MAKER_EXE" "$MAKER_OPTS" "$MAKER_BOPTS" \
    -qq -base "$MAKERDIR" -fix_nucleotides

check_maker_output "Round 2"

gff3_merge \
    -d "${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log" \
    -o "round2/${genotype}.snap02.all.gff"

log "  Round 2 gene models: $(grep -c $'\tgene\t' "round2/${genotype}.snap02.all.gff" || true)"

# SNAP Training — Round 2
log "=== SNAP training on Round 2 predictions ==="
mkdir -p snap.round2 train_snap.round2

(
    cd snap.round2
    maker2zff "../round2/${genotype}.snap02.all.gff"
    [[ -s genome.ann && -s genome.dna ]] || { echo "ERROR: maker2zff did not produce genome.ann/genome.dna" >&2; exit 1; }
    fathom genome.ann genome.dna -validate > snap_validate_output.txt
    grep 'errors' snap_validate_output.txt | awk '{print $2}' > snap.gene_model_error.txt || true
    grep -vwE -f snap.gene_model_error.txt genome.ann > genome.ann1 || cp genome.ann genome.ann1
)

(
    cd train_snap.round2
    fathom ../snap.round2/genome.ann1 ../snap.round2/genome.dna -categorize 1000
    fathom uni.ann uni.dna -export 1000 -plus
    forge export.ann export.dna
    hmm-assembler.pl "${genotype}.snap2" . > "${genotype}.snap2.hmm"
)

[[ -f "train_snap.round2/${genotype}.snap2.hmm" ]] || err "SNAP Round 2 HMM not produced."
log "  SNAP Round 2 HMM: train_snap.round2/${genotype}.snap2.hmm"

# ===========================================================================
# ROUND 3 — SNAP (Round 2 HMM)
# maker_opts.ctl: maker_gff=<abs>/round2/<genotype>.snap02.all.gff,
#   snaphmm=<abs>/train_snap.round2/<genotype>.snap2.hmm
# ===========================================================================
log "=== Round 3: SNAP Round 2 HMM ==="
mkdir -p round3

mpiexec -n "$cpus" maker "$MAKER_EXE" "$MAKER_OPTS" "$MAKER_BOPTS" \
    -qq -base "$MAKERDIR" -fix_nucleotides

check_maker_output "Round 3"

gff3_merge \
    -d "${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log" \
    -o "round3/${genotype}.snap03.all.gff"

log "  Round 3 gene models: $(grep -c $'\tgene\t' "round3/${genotype}.snap03.all.gff" || true)"

# SNAP Training — Round 3 (3 rounds total for convergence)
log "=== SNAP training on Round 3 predictions ==="
mkdir -p snap.round3 train_snap.round3

(
    cd snap.round3
    maker2zff "../round3/${genotype}.snap03.all.gff"
    [[ -s genome.ann && -s genome.dna ]] || { echo "ERROR: maker2zff did not produce genome.ann/genome.dna" >&2; exit 1; }
    fathom genome.ann genome.dna -validate > snap_validate_output.txt
    grep 'errors' snap_validate_output.txt | awk '{print $2}' > snap.gene_model_error.txt || true
    grep -vwE -f snap.gene_model_error.txt genome.ann > genome.ann1 || cp genome.ann genome.ann1
)

(
    cd train_snap.round3
    fathom ../snap.round3/genome.ann1 ../snap.round3/genome.dna -categorize 1000
    fathom uni.ann uni.dna -export 1000 -plus
    forge export.ann export.dna
    hmm-assembler.pl "${genotype}.snap3" . > "${genotype}.snap3.hmm"
)

[[ -f "train_snap.round3/${genotype}.snap3.hmm" ]] || err "SNAP Round 3 HMM not produced."
log "  SNAP Round 3 HMM: train_snap.round3/${genotype}.snap3.hmm"

# ===========================================================================
# ROUND 4 (FINAL) — SNAP Round 3 + GeneMark (BRAKER2) + Augustus
# maker_opts.ctl: maker_gff=<abs>/round3/<genotype>.snap03.all.gff,
#   snaphmm=<abs>/train_snap.round3/<genotype>.snap3.hmm,
#   gmhmm=<abs>/braker/GeneMark_hmm_full.mod, augustus_species=ensete_ventricosum,
#   est2genome=0, protein2genome=0
# ===========================================================================
log "=== Round 4 (final): SNAP Round 3 + GeneMark + Augustus ==="
mkdir -p round4_final

mpiexec -n "$cpus" maker "$MAKER_EXE" "$MAKER_OPTS" "$MAKER_BOPTS" \
    -qq -base "$MAKERDIR" -fix_nucleotides

check_maker_output "Round 4 (final)"

# Merge GFF3 and FASTA across all scaffolds
FINAL_GFF="round4_final/${genotype}.final.all.gff"
FINAL_PROTEINS="round4_final/${genotype}.final.proteins.fasta"
FINAL_TRANSCRIPTS="round4_final/${genotype}.final.transcripts.fasta"

gff3_merge \
    -d "${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log" \
    -o "$FINAL_GFF"

fasta_merge \
    -d "${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log" \
    -o "round4_final/${genotype}.final"

[[ -f "$FINAL_GFF" ]] || err "Final merged GFF not produced: $FINAL_GFF"

FINAL_GENE_COUNT=$(grep -c $'\tgene\t' "$FINAL_GFF" || true)
log "  Final gene models: $FINAL_GENE_COUNT"

# AED distribution QC (AED=0 perfect agreement; high-quality annotation >80% genes with AED<0.5)
if [[ -n "$AED_CDF" && -x "$AED_CDF" ]]; then
    "$AED_CDF" -b 0.01 "$FINAL_GFF" 2>/dev/null | head -20 | \
        while IFS= read -r line; do log "    $line"; done
else
    log "  WARNING: AED_cdf_generator.pl not found — run manually: AED_cdf_generator.pl -b 0.01 $FINAL_GFF"
fi

# Functional annotation — uncomment and configure once BLAST/InterPro results are ready:
#   maker_functional_fasta -g <genome.fasta> -f $FINAL_PROTEINS -b <blast_output.txt> \
#       > round4_final/${genotype}.final.proteins.functional.fasta
#   maker_functional_gff  -g <genome.fasta> -f $FINAL_GFF        -b <blast_output.txt> \
#       > round4_final/${genotype}.final.functional.gff
log "NOTE: Run maker_functional_fasta / maker_functional_gff after BLAST vs UniProt/SwissProt."

log "=== MAKER iterative annotation pipeline complete ==="
log "  Final GFF3:      $FINAL_GFF ($FINAL_GENE_COUNT genes)"
log "  Proteins FASTA:  round4_final/${genotype}.final.all.maker.proteins.fasta"
log "  Transcripts:     round4_final/${genotype}.final.all.maker.transcripts.fasta"
log "  SNAP HMM (R3):   train_snap.round3/${genotype}.snap3.hmm"
