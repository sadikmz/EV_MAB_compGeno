#!/bin/bash
# =============================================================================
# visualize_alphafold_results.sh
# Post-processing and visualisation of AlphaFold2 structural predictions for
# EV/MAB novel (orphan) proteins with no detectable homologues.
#
# The Python helper script produces:
#   - Per-model pLDDT score plots: per-residue Local Distance Difference Test
#     confidence (0–100); >90 = very high, 70–90 = high, 50–70 = low,
#     <50 = very low (disordered/unreliable).
#   - Predicted Aligned Error (PAE) matrices: 2-D plots showing inter-residue
#     distance error; low PAE (dark blue) indicates confident relative domain
#     positioning.
#   - Structure confidence summaries: mean pLDDT, % residues > 70, PAE range.
#   - Model ranking: ranked_0.pdb … ranked_4.pdb by model confidence score.
# =============================================================================

set -euo pipefail

# ===========================================================================
# Helper functions
# ===========================================================================

log() {
    # Structured timestamped info messages — captured in SLURM .out file.
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"
}

err() {
    # Structured timestamped error messages — captured in SLURM .err file.
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

check_tool() {
    # Verify that a required executable is accessible before the main run.
    local tool="$1"
    if ! command -v "${tool}" &>/dev/null; then
        err "Required tool not found on PATH: ${tool}"
        err "Ensure the appropriate module is loaded."
        exit 1
    fi
}

# ===========================================================================
# Named variable block
# ===========================================================================

# --- Protein/query name: accept as first positional argument or use default.
#     Usage: bash visualize_alphafold_results.sh [NAME]
#     NAME should match the subdirectory created by AlphaFold for this protein.
NAME="${1:-rbd}"

# --- AlphaFold output directory layout: <NAME>/<NAME>/
#     AlphaFold writes per-query results to <output_dir>/<query_name>/
DIR="${NAME}/${NAME}/"

# --- Python visualisation script (must reside alongside this shell script)
VIZ_SCRIPT="$(dirname "${BASH_SOURCE[0]}")/visualize_alphafold_results.py"

# ===========================================================================
# Module environment
# ===========================================================================

module purge
# Python 3.10 with GCC 12.2.0 toolchain — matches AlphaFold conda/module env.
module load Python/3.10.8-GCCcore-12.2.0
# matplotlib: required for pLDDT line plots and PAE heatmaps.
module load matplotlib

# ===========================================================================
# Tool availability checks
# ===========================================================================

check_tool python

# Capture Python version for Methods reproducibility.
log "=== Tool versions ==="
log "Python: $(python --version 2>&1)"
log "matplotlib: $(python -c 'import matplotlib; print(matplotlib.__version__)' 2>&1)"

# ===========================================================================
# Input validation
# ===========================================================================

log "=== Input validation ==="
log "Query name      : ${NAME}"
log "AlphaFold dir   : ${DIR}"

# Verify that the AlphaFold results directory exists.
if [[ ! -d "${DIR}" ]]; then
    err "AlphaFold results directory not found: ${DIR}"
    err "Expected layout: <NAME>/<NAME>/ — check that AlphaFold completed for '${NAME}'."
    exit 1
fi

# AlphaFold writes 5 ranked PDB models (ranked_0.pdb … ranked_4.pdb).
# Require at least one .pdb before attempting visualisation.
n_pdb=$(find "${DIR}" -maxdepth 1 -name "*.pdb" | wc -l)
if [[ "${n_pdb}" -eq 0 ]]; then
    err "No .pdb files found in ${DIR}"
    err "AlphaFold may not have completed successfully for '${NAME}'."
    exit 1
fi

log "Found ${n_pdb} PDB model file(s) in ${DIR}"

# Verify the presence of AlphaFold confidence JSON files (contain pLDDT/PAE).
n_json=$(find "${DIR}" -maxdepth 1 -name "*.json" | wc -l)
log "Found ${n_json} JSON confidence file(s) in ${DIR}"

# Ensure the Python visualisation script itself is present.
if [[ ! -f "${VIZ_SCRIPT}" ]]; then
    err "Visualisation script not found: ${VIZ_SCRIPT}"
    err "Place visualize_alphafold_results.py in the same directory as this script."
    exit 1
fi

log "Visualisation script: ${VIZ_SCRIPT}"

# ===========================================================================
# Run AlphaFold visualisation
# ===========================================================================

log "=== Starting AlphaFold visualisation for '${NAME}' ==="
log "Input  directory: ${DIR}"
log "Output directory: ${DIR}  (figures written alongside predictions)"

python "${VIZ_SCRIPT}" \
    --input_dir "${DIR}" \
    --output_dir "${DIR}"
    # Note: --name flag is available in the script but omitted here because
    # the input/output directory already scopes output to this query protein.
    # Uncomment and extend if the script supports per-model naming:
    # --name "${NAME}"

# ===========================================================================
# Post-run summary
# ===========================================================================

log "=== Post-run summary for '${NAME}' ==="

# Count output figures (PNG/PDF) written by the visualisation script.
n_png=$(find "${DIR}" -maxdepth 1 -name "*.png" | wc -l)
n_pdf=$(find "${DIR}" -maxdepth 1 -name "*.pdf" | wc -l)
log "PDB model files : ${n_pdb}"
log "Output PNG plots: ${n_png}"
log "Output PDF plots: ${n_pdf}"

# List every output file for the SLURM log — useful for provenance tracking.
log "Output files produced:"
find "${DIR}" -maxdepth 1 \( -name "*.png" -o -name "*.pdf" \) \
    | sort | while IFS= read -r f; do
        log "  ${f}"
    done

# Report the top-ranked model file for downstream structural analysis.
top_model="${DIR}/ranked_0.pdb"
if [[ -f "${top_model}" ]]; then
    log "Top-ranked model (highest confidence): ${top_model}"
else
    log "ranked_0.pdb not found; check AlphaFold output integrity."
fi

log "=== AlphaFold visualisation complete for '${NAME}' ==="
