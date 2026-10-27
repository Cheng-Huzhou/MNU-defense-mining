#!/usr/bin/env bash
set -euo pipefail

HITS_FASTA=${1:?usage: 11_cluster_align_tree.sh mu_hits.faa [prefix] [threads] [model] [iqtree_flags]}
PREFIX=${2:-mu}
THREADS=${3:-32}
MODEL=${4:-MFP}
IQFLAGS=${5:-}

mmseqs easy-cluster "$HITS_FASTA" "${PREFIX}_clu" "${PREFIX}_tmp" \
    --min-seq-id 0.5 \
    -c 0.7 \
    --cov-mode 1 \
    --threads "$THREADS"

mafft --auto --thread "$THREADS" "${PREFIX}_clu_rep_seq.fasta" > "${PREFIX}_aln.fasta"
iqtree -s "${PREFIX}_aln.fasta" -m "$MODEL" -T AUTO $IQFLAGS --prefix "${PREFIX}_tree"
