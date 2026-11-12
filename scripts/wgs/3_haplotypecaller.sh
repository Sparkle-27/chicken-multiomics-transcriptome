#!/usr/bin/env bash
# Per-sample variant calling with GATK HaplotypeCaller.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 GRCg7b.fa input.bam output.vcf.gz" >&2
    exit 2
fi
reference=$1
bam=$2
output=$3
mkdir -p "$(dirname "$output")"
gatk --java-options -Xmx4G HaplotypeCaller -R "$reference" -I "$bam" -O "$output"
