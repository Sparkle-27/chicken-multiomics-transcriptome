#!/usr/bin/env bash
# Provide the matching bulk FLNC input.
set -euo pipefail
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 clustered.bam GRCg7b.fa flnc_input_for_collapse threads output_prefix" >&2
    exit 2
fi
prefix=$5
isoseq_bin=${ISOSEQ_BIN:-isoseq3}
mkdir -p "$(dirname "$prefix")"
pbmm2 align --preset ISOSEQ --sort "$1" "$2" "${prefix}-pbmm2.bam" -j "$4"
"$isoseq_bin" collapse --do-not-collapse-extra-5exons --max-5p-diff 200 \
    "${prefix}-pbmm2.bam" "$3" "${prefix}.gff" -j "$4"
