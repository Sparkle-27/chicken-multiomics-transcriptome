#!/usr/bin/env bash
# Reference-guided StringTie quantification with -e.
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 input.bam Ensembl110.gtf threads output_prefix" >&2
    exit 2
fi
prefix=$4
mkdir -p "$(dirname "$prefix")"
stringtie -e "$1" -p "$3" -G "$2" -o "${prefix}.gtf" -A "${prefix}.tab"
