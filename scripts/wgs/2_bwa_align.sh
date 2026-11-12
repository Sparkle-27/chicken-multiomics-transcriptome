#!/usr/bin/env bash
# BWA-MEM alignment, coordinate sorting and mapping statistics.
set -euo pipefail
if [[ $# -ne 6 ]]; then
    echo "Usage: $0 GRCg7b.fa clean_R1.fq.gz clean_R2.fq.gz sample_id threads output_prefix" >&2
    exit 2
fi
reference=$1
r1=$2
r2=$3
sample=$4
threads=$5
prefix=$6
mkdir -p "$(dirname "$prefix")"
rg="@RG\tID:${sample}\tLB:${sample}\tPL:ILLUMINA\tSM:${sample}"
bwa mem -t "$threads" -M -R "$rg" "$reference" "$r1" "$r2" |
    samtools sort -@ "$threads" -o "${prefix}.sort.bam" -
samtools index "${prefix}.sort.bam"
samtools flagstat "${prefix}.sort.bam" > "${prefix}.flagstat.txt"
