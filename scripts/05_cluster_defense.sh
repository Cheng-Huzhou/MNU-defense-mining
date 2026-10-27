#!/usr/bin/env bash
set -euo pipefail

PROTEINS=${1:?usage: 05_cluster_defense.sh proteins.faa [out_prefix] [threads]}
PREFIX=${2:-defense_clusters}
THREADS=${3:-32}

mmseqs easy-cluster "$PROTEINS" "$PREFIX" "${PREFIX}_tmp" \
    --min-seq-id 0.8 \
    -c 0.8 \
    --cov-mode 1 \
    --threads "$THREADS"
