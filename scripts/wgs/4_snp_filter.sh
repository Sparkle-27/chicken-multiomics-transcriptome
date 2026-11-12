#!/usr/bin/env bash
# Annotate SNP FILTER fields using QD, FS, SOR, MQ and rank-sum thresholds.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 GRCg7b.fa SNP_only_input.vcf.gz SNP_tagged_output.vcf.gz" >&2
    exit 2
fi
mkdir -p "$(dirname "$3")"
gatk VariantFiltration -R "$1" -V "$2" -O "$3" \
    --filter-name QD_lt_5 --filter-expression 'QD < 5.0' \
    --filter-name FS_gt_60 --filter-expression 'FS > 60.0' \
    --filter-name SOR_gt_3 --filter-expression 'SOR > 3.0' \
    --filter-name MQ_lt_40 --filter-expression 'MQ < 40.0' \
    --filter-name MQRankSum_lt_neg12_5 --filter-expression 'MQRankSum < -12.5' \
    --filter-name ReadPosRankSum_lt_neg8 --filter-expression 'ReadPosRankSum < -8.0'
