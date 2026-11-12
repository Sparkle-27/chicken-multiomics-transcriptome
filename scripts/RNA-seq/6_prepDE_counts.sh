#!/usr/bin/env bash
# Gene and transcript counts from StringTie GTF files.
set -euo pipefail
if [[ $# -ne 4 || ! ${3:-} =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 /path/to/prepDE.py3 sample_gtf_list.txt average_read_length output_prefix" >&2
    exit 2
fi
prefix=$4
mkdir -p "$(dirname "$prefix")"
python3 "$1" -i "$2" -l "$3" -g "${prefix}.gene_counts.csv" -t "${prefix}.transcript_counts.csv"
