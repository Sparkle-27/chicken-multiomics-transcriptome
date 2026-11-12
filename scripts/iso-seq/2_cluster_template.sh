#!/usr/bin/env bash
# Methods and the commented cluster2 command in script/21-ccs2Isoform-s1.sh.
set -euo pipefail
if [[ $# -ne 4 || ( ${1:-} != cluster && ${1:-} != cluster2 ) ]]; then
    echo "Usage: $0 confirmed_subcommand(cluster|cluster2) flnc_input threads output.clustered.bam" >&2
    exit 2
fi
isoseq_bin=${ISOSEQ_BIN:-isoseq3}
mkdir -p "$(dirname "$4")"
"$isoseq_bin" "$1" "$2" "$4" -j "$3"
