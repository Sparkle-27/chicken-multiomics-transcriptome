#!/usr/bin/env bash
# Mark duplicates after name sorting and fixmate; retain marked reads.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 input.bam threads output.markdup.bam" >&2
    exit 2
fi
bam=$1
threads=$2
output=$3
mkdir -p "$(dirname "$output")"
samtools sort -n -u -@ "$threads" "$bam" |
    samtools fixmate -m -u - - |
    samtools sort -u -@ "$threads" - |
    samtools markdup -s -f "${output}.stats.txt" - "$output"
samtools index "$output"
