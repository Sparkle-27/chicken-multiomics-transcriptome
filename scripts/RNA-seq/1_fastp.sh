#!/usr/bin/env bash
# Paired-end read filtering with fastp 0.23.2 defaults.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 R1.fq.gz R2.fq.gz output_prefix" >&2
    exit 2
fi
r1=$1
r2=$2
prefix=$3
mkdir -p "$(dirname "$prefix")"
fastp -i "$r1" -I "$r2" \
    -o "${prefix}_R1.fq.gz" -O "${prefix}_R2.fq.gz" \
    --html "${prefix}.fastp.html" --json "${prefix}.fastp.json"
