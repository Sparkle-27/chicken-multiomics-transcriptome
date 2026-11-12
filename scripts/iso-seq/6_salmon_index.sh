#!/usr/bin/env bash
# Adapted from script/22-isoformQuant.sh. Use the confirmed high-confidence GTF.
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 final_catalog.gtf GRCg7b.fa threads output_prefix" >&2
    exit 2
fi
prefix=$4
mkdir -p "$(dirname "$prefix")"
gffread "$1" -g "$2" -w "${prefix}.fa"
salmon index -t "${prefix}.fa" -p "$3" -i "${prefix}.salmon_index"
