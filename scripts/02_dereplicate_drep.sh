#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR=${1:?usage: 02_dereplicate_drep.sh genome_dir [outdir] [threads]}
OUTDIR=${2:-drep_out}
THREADS=${3:-32}

dRep dereplicate "$OUTDIR" \
    -g "$GENOME_DIR"/*.fna \
    -p "$THREADS" \
    -pa 0.9 \
    -sa 0.99 \
    -nc 0.1 \
    --ignoreGenomeQuality

echo "representatives: $OUTDIR/dereplicated_genomes"
