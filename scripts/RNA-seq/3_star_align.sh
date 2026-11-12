#!/usr/bin/env bash
# Paired-end STAR alignment and sorted BAM output.
set -euo pipefail
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 STAR_index clean_R1.fq.gz clean_R2.fq.gz threads output_prefix" >&2
    exit 2
fi
index=$1
r1=$2
r2=$3
threads=$4
prefix=$5
mkdir -p "$(dirname "$prefix")"
STAR --runMode alignReads --runThreadN "$threads" --genomeDir "$index" \
    --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
    --outSJfilterOverhangMin 30 16 16 16 --outSJfilterCountUniqueMin 4 2 2 2 \
    --alignSJoverhangMin 6 --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --alignIntronMax 500000 --alignMatesGapMax 500000 \
    --outFilterMismatchNoverReadLmax 0.07 \
    --outFileNamePrefix "${prefix}." --outStd SAM --outSAMmode Full |
    samtools sort -@ "$threads" -o "${prefix}.sort.bam" -
samtools index "${prefix}.sort.bam"
