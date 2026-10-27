#!/usr/bin/env python3
import argparse
import collections
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-r", "--adjacency", required=True)
    ap.add_argument("-t", "--taxonomy", required=True)
    ap.add_argument("-o", "--output", default="phylum_distribution.tsv")
    args = ap.parse_args()

    genomes = set()
    with open(args.adjacency, newline="") as f:
        for row in csv.DictReader(f):
            genomes.add(row["genome"])

    phylum_by_genome = {}
    with open(args.taxonomy, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            phylum_by_genome[row["genome"]] = row.get("phylum", "") or "unassigned"

    counts = collections.Counter(phylum_by_genome.get(g, "unassigned") for g in genomes)
    total = sum(counts.values()) or 1

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["phylum", "count", "percent"])
        for phylum, n in counts.most_common():
            w.writerow([phylum, n, f"{100.0 * n / total:.1f}"])

    print(f"genomes: {total}, phyla: {len(counts)}")


if __name__ == "__main__":
    main()
