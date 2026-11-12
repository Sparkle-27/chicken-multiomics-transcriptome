#!/usr/bin/env bash
# Adapted from 11-script/22-SQANTI3-QC-s1-v2.bash. Supply the actual study rules JSON.
set -euo pipefail
if [[ $# -ne 9 ]]; then
    echo "Usage: $0 SQANTI3_directory isoforms.gtf Ensembl110.gtf GRCg7b.fa FL_counts short_reads_list rules.json threads output_prefix" >&2
    exit 2
fi
sqanti=$1
isoforms=$2
reference_gtf=$3
reference_fasta=$4
fl_counts=$5
short_reads=$6
rules=$7
threads=$8
prefix=$9
if [[ ! -s "$rules" ]]; then
    echo "Provide the confirmed, nonempty SQANTI3 rules JSON." >&2
    exit 2
fi
out_dir=$(dirname "$prefix")
name=$(basename "$prefix")
mkdir -p "$out_dir"
sqanti=$(cd "$sqanti" && pwd)
rules="$(cd "$(dirname "$rules")" && pwd)/$(basename "$rules")"
python3 "$sqanti/sqanti3_qc.py" "$isoforms" "$reference_gtf" "$reference_fasta" \
    -o "$name" -d "$out_dir" -fl "$fl_counts" --short_reads "$short_reads" \
    --isoform_hits --saturation --cpus "$threads" --report both
cd "$out_dir"
python3 "$sqanti/sqanti3_filter.py" rules \
    --isoforms "${name}_corrected.fasta" --gtf "${name}_corrected.gtf" \
    --faa "${name}_corrected.faa" --json_filter "$rules" "${name}_classification.txt"
