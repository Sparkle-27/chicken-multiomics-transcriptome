#!/usr/bin/env bash
# Adapted from script/21-ccs2Isoform-s1.sh. This runs CCS on the complete input BAM.
set -euo pipefail
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 subreads.bam primers.fasta primer_pair_label threads output_prefix" >&2
    exit 2
fi
subreads=$1
primers=$2
primer_pair=$3
threads=$4
prefix=$5
isoseq_bin=${ISOSEQ_BIN:-isoseq3}
mkdir -p "$(dirname "$prefix")"
ccs --hifi-kinetics --min-passes 3 --min-rq 0.99 --max-length 100000 \
    --min-length 100 "$subreads" "${prefix}.ccs.bam" -j "$threads"
lima "${prefix}.ccs.bam" "$primers" "${prefix}.fl.bam" --isoseq --peek-guess -j "$threads"
"$isoseq_bin" refine "${prefix}.fl.${primer_pair}.bam" "$primers" \
    "${prefix}.flnc.bam" --require-polya
