#!/usr/bin/env bash
# The v2 source records a full rules rescue. Confirm this branch and the study JSON.
set -euo pipefail
if [[ $# -ne 6 ]]; then
    echo "Usage: $0 SQANTI3_directory qc_output_prefix Ensembl110.gtf GRCg7b.fa rules.json rescue_output_prefix" >&2
    exit 2
fi
prefix=$2
rules=$5
output=$6
if [[ ! -s "$rules" ]]; then
    echo "Provide the confirmed, nonempty SQANTI3 rules JSON." >&2
    exit 2
fi
mkdir -p "$(dirname "$output")"
python3 "$1/sqanti3_rescue.py" rules \
    --isoforms "${prefix}_corrected.fasta" --gtf "${prefix}.filtered.gtf" \
    -g "$3" -f "$4" -k "${prefix}_classification.txt" -e all \
    --mode full -j "$rules" -o "$(basename "$output")" -d "$(dirname "$output")" \
    "${prefix}_RulesFilter_result_classification.txt"
