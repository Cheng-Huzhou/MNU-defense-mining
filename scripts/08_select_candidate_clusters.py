#!/usr/bin/env python3
import argparse
import collections
import csv
import glob
import os


def load_known_systems(data_dir):
    known = collections.defaultdict(set)

    for sys_tsv in glob.glob(os.path.join(data_dir, "*", "defensefinder", "*.tsv")):
        if not sys_tsv.endswith("_systems.tsv"):
            continue
        with open(sys_tsv, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                system = row.get("type", "")
                for pid in row.get("protein_in_syst", "").split(","):
                    if pid:
                        known[pid].add(f"df:{system}")

    for padloc_csv in glob.glob(os.path.join(data_dir, "*", "padloc", "*_padloc.csv")):
        genome_dir = os.path.dirname(os.path.dirname(padloc_csv))
        coord_to_protein = {}
        for defense_csv in glob.glob(os.path.join(genome_dir, "*_defense.csv")):
            with open(defense_csv, newline="") as f:
                for row in csv.DictReader(f):
                    key = (row.get("genomic_accession"), row.get("start"), row.get("end"), row.get("strand"))
                    pid = row.get("product_accession")
                    if pid:
                        coord_to_protein[key] = pid
        with open(padloc_csv, newline="") as f:
            for row in csv.DictReader(f):
                key = (row.get("seqid"), row.get("start"), row.get("end"), row.get("strand"))
                pid = coord_to_protein.get(key)
                if pid:
                    known[pid].add(f"pd:{row.get('system', '')}")

    return known


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--clusters", required=True)
    ap.add_argument("-d", "--data-dir", required=True)
    ap.add_argument("-e", "--ecoli-ids", default="ecoli_hit_ids.txt")
    ap.add_argument("-o", "--output", default="candidate_clusters.tsv")
    args = ap.parse_args()

    known = load_known_systems(args.data_dir)
    ecoli_ids = set()
    if os.path.exists(args.ecoli_ids):
        ecoli_ids = {line.strip() for line in open(args.ecoli_ids) if line.strip()}

    clusters = collections.OrderedDict()
    with open(args.clusters, newline="") as f:
        for line in f:
            if not line.strip():
                continue
            rep, member = line.rstrip("\n").split("\t")[:2]
            clusters.setdefault(rep, []).append(member)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rep_protein", "cluster_size", "status", "known_systems",
                    "ecoli_homolog", "category"])
        for rep, members in clusters.items():
            systems = known.get(rep, set())
            status = "known_defense" if systems else "novel"
            has_ecoli = "yes" if rep in ecoli_ids else "no"
            if status == "known_defense":
                category = "known_defense"
            elif has_ecoli == "yes":
                category = "novel_ecoli_conserved"
            else:
                category = "novel_no_ecoli_homolog"
            w.writerow([rep, len(members), status, ";".join(sorted(systems)),
                        has_ecoli, category])


if __name__ == "__main__":
    main()
