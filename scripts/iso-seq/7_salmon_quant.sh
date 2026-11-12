#!/usr/bin/env bash
# Adapted from script/22-isoformQuant.sh. Use the catalog for the matching tissue.
set -euo pipefail
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 salmon_index clean_R1.fq.gz clean_R2.fq.gz threads output_dir" >&2
    exit 2
fi
mkdir -p "$(dirname "$5")"
salmon quant -i "$1" --libType A -1 "$2" -2 "$3" \
    --validateMappings -p "$4" -o "$5"
