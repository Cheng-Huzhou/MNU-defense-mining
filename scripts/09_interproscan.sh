#!/usr/bin/env bash
set -euo pipefail

PROTEINS=${1:?usage: 09_interproscan.sh proteins.faa [output.tsv] [cpu] [appl]}
OUTPUT=${2:-interproscan.tsv}
CPU=${3:-16}
APPL=${4:-Pfam}

interproscan.sh -i "$PROTEINS" -f tsv -dp -appl "$APPL" -o "$OUTPUT" -cpu "$CPU"
