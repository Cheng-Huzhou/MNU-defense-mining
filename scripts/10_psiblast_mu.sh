#!/usr/bin/env bash
set -euo pipefail

MU=${1:?usage: 10_psiblast_mu.sh MU.faa nr_db [threads]}
NR=${2:?}
THREADS=${3:-32}

psiblast -query "$MU" -db "$NR" \
    -evalue 1e-5 \
    -num_iterations 3 \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen" \
    -out mu_hits.tsv
