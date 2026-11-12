#!/usr/bin/env bash
# STAR genome index: sjdbOverhang 149, genomeSAindexNbases 13.
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 GRCg7b.fa Ensembl110.gtf index_directory threads" >&2
    exit 2
fi
mkdir -p "$3"
STAR --runMode genomeGenerate --runThreadN "$4" --genomeDir "$3" \
    --genomeFastaFiles "$1" --sjdbGTFfile "$2" \
    --genomeSAindexNbases 13 --sjdbOverhang 149
