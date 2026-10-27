#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=${1:?usage: 06_annotate_known.sh ncbi_dataset_data_dir [threads]}
THREADS=${2:-16}

for d in "$DATA_DIR"/*/; do
    faa=$(find "$d" -name "*_protein.faa" | head -1)
    fna=$(find "$d" -name "*_genomic.fna" ! -name "*_cds_from_genomic.fna" | head -1)

    if [ -n "$faa" ]; then
        conda run -n defensefinder defense-finder run "$faa" -o "$d/defensefinder" -w "$THREADS"
    fi
    if [ -n "$fna" ]; then
        mkdir -p "$d/padloc"
        conda run -n padloc padloc --fna "$fna" --outdir "$d/padloc"
    fi
done
