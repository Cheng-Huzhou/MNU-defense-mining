#!/usr/bin/env bash
set -euo pipefail

ACCESSION_FILE=${1:?usage: 01_download_assemblies.sh accessions.txt [outdir]}
OUTDIR=${2:-genomes}

mkdir -p "$OUTDIR"
datasets download genome accession \
    --inputfile "$ACCESSION_FILE" \
    --include genome,protein,gff3,cds \
    --filename "$OUTDIR/assemblies.zip"
unzip -q -o "$OUTDIR/assemblies.zip" -d "$OUTDIR/assemblies"

echo "genomes: $OUTDIR/assemblies/ncbi_dataset/data"
