#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------
# paf2bed.sh: Convert a single PAF file to merged BEDs for query and reference
# Usage: ./paf2bed.sh INPUT.paf QUERY_GENOME_NAME REFERENCE_GENOME_NAME
# Outputs: query_aligned.bed and reference_aligned.bed
# -----------------------------------------------

# Check for required arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 INPUT.paf QUERY_GENOME_NAME REFERENCE_GENOME_NAME" >&2
    exit 1
fi

input_paf="$1"   # User-provided PAF file
query="$2"       # User-provided query genome name
ref="$3"         # User-provided reference genome name

# Check if input file exists
if [[ ! -f "$input_paf" ]]; then
    echo "Error: Input file '$input_paf' not found!" >&2
    exit 1
fi

# Check if bedtools is installed and available in PATH
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is required but not found in PATH" >&2
    exit 1
fi

# Temporary files for intermediate BEDs
tmp_query_bed=$(mktemp)
tmp_ref_bed=$(mktemp)

# Clean up temp files on exit
trap 'rm -f "$tmp_query_bed" "$tmp_ref_bed"' EXIT

# Extract query BED fields: columns 1 (query), 3 (q_start), 4 (q_end), 10, 11
cut -f1,3,4,10,11 "$input_paf" > "$tmp_query_bed"

# Extract reference BED fields: columns 6 (ref), 8 (r_start), 9 (r_end), 10, 11
cut -f6,8,9,10,11 "$input_paf" > "$tmp_ref_bed"

# Function to process a BED file:
# - Ensures start < end, sorts, merges, and sums columns 4 and 5
# - Outputs to the specified output file
process_bed() {
    local input=$1      # Input BED file
    local output=$2     # Output BED file

    awk -v OFS='\t' '
        {
            if ($2 > $3) print $1, $3, $2, $4, $5;
            else print $1, $2, $3, $4, $5
        }' "$input" |
    bedtools sort |
    bedtools merge |
    awk -v OFS='\t' '{print $1, $2, $3, $3-$2 }' > "$output"
}

# Process and merge for query
process_bed "$tmp_query_bed" "$2"_aligned.bed

# Process and merge for reference
process_bed "$tmp_ref_bed" "$3"_aligned.bed

