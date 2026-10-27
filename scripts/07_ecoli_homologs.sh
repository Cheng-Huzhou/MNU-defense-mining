#!/usr/bin/env bash
set -euo pipefail

QUERY=${1:?usage: 07_ecoli_homologs.sh candidates.faa ecoli_proteins.faa [threads]}
ECOLI_FASTA=${2:?}
THREADS=${3:-32}

mmseqs createdb "$ECOLI_FASTA" ecoli_db
mmseqs easy-search "$QUERY" ecoli_db ecoli_hits.m8 ecoli_tmp \
    -e 1e-5 \
    -s 7 \
    --threads "$THREADS"

cut -f1 ecoli_hits.m8 | sort -u > ecoli_hit_ids.txt
echo "queries with E. coli homolog: $(wc -l < ecoli_hit_ids.txt)"
