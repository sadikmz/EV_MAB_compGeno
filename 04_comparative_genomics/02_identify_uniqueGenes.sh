#!/bin/bash

set -euo pipefail

# This analysis is implemented as a Python pipeline.
# See: 02_identify_uniqueGenes/pav_analysis.py
#
# Usage example:
#   python 02_identify_uniqueGenes/pav_analysis.py \
#     --query-genome EVMZ.softmasked.fna \
#     --query-prefix EVMZ \
#     --query-gff EVMZ.genes.gff3 \
#     --ref-genome MA.fna \
#     --ref-prefix MA \
#     --ref-gff MA.genes.gff3 \
#     --sra-ids SRR27722736 \
#     --output-dir pav_output/ \
#     --cpus 28
