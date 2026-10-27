#!/usr/bin/env python3
import argparse
import csv


def base_accession(name):
    parts = name.split("_")
    return f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--accessions", required=True)
    ap.add_argument("-t", "--gtdb-taxonomy", required=True)
    ap.add_argument("-o", "--output", default="gtdb_phylum.tsv")
    args = ap.parse_args()

    phylum_by_accession = {}
    with open(args.gtdb_taxonomy) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            leaf = fields[0]
            if leaf.startswith("RS_") or leaf.startswith("GB_"):
                leaf = leaf[3:]
            ranks = dict(x.split("__", 1) for x in fields[1].split(";") if "__" in x)
            phylum_by_accession[leaf] = ranks.get("p", "")

    with open(args.accessions) as f:
        accessions = [line.strip() for line in f if line.strip()]

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["genome", "phylum"])
        for acc in accessions:
            w.writerow([acc, phylum_by_accession.get(base_accession(acc), "")])


if __name__ == "__main__":
    main()
