#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=${1:?usage: 03_predict_defense.sh ncbi_dataset_data_dir}

for d in "$DATA_DIR"/*/; do
    acc=$(basename "$d")
    ft=$(find "$d" -name "*_feature_table.txt" | head -1)
    cds=$(find "$d" -name "*_cds_from_genomic.fna" | head -1)
    faa=$(find "$d" -name "*_protein.faa" | head -1)

    if [ -z "$faa" ]; then
        echo "skip $acc: no protein.faa" >&2
        continue
    fi

    defense_predictor \
        --ncbi_feature_table "$ft" \
        --ncbi_cds_from_genomic "$cds" \
        --ncbi_protein_fasta "$faa" \
        --output "$d/${acc}_defense.csv"
done
