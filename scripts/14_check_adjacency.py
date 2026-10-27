#!/usr/bin/env python3
import argparse
import csv
import gzip
import os


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_gff_coords(gff_path, target_ids):
    coords = {}
    with open_text(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            protein_id = None
            for attr in parts[8].split(";"):
                if attr.startswith("protein_id="):
                    protein_id = attr.split("=", 1)[1]
                    break
            if protein_id in target_ids:
                coords[protein_id] = {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
    return coords


def distance(mu, nu):
    if mu["contig"] != nu["contig"]:
        return None
    if mu["end"] < nu["start"]:
        return nu["start"] - mu["end"]
    if nu["end"] < mu["start"]:
        return mu["start"] - nu["end"]
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-r", "--hits", required=True)
    ap.add_argument("-g", "--genome-dir", required=True)
    ap.add_argument("-m", "--mu-id", required=True)
    ap.add_argument("-d", "--max-distance", type=int, default=10000)
    ap.add_argument("-o", "--output", default="adjacent_MU_NU.csv")
    args = ap.parse_args()

    hits = {}
    with open(args.hits, newline="") as f:
        for row in csv.DictReader(f):
            g = row["genome"]
            hits.setdefault(g, {"MU": set(), "NU": set()})
            hits[g]["MU" if row["qseqid"] == args.mu_id else "NU"].add(row["subject"])

    pairs = []
    for genome, ids in hits.items():
        if not ids["MU"] or not ids["NU"]:
            continue
        gff = None
        for suffix in ("_genomic.gff.gz", "_genomic.gff", ".gff"):
            path = os.path.join(args.genome_dir, genome, f"{genome}{suffix}")
            if os.path.exists(path):
                gff = path
                break
        if gff is None:
            gffs = [p for p in os.listdir(os.path.join(args.genome_dir, genome)) if p.endswith(".gff") or p.endswith(".gff.gz")]
            if not gffs:
                continue
            gff = os.path.join(args.genome_dir, genome, gffs[0])

        coords = parse_gff_coords(gff, ids["MU"] | ids["NU"])
        mu_coords = {k: coords[k] for k in ids["MU"] if k in coords}
        nu_coords = {k: coords[k] for k in ids["NU"] if k in coords}

        for mu_id, mu in mu_coords.items():
            for nu_id, nu in nu_coords.items():
                dist = distance(mu, nu)
                if dist is not None and dist <= args.max_distance:
                    pairs.append([genome, mu_id, nu_id, mu["contig"], dist,
                                  mu["strand"] == nu["strand"]])

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome", "mu_id", "nu_id", "contig", "distance", "same_strand"])
        w.writerows(pairs)

    print(f"adjacent pairs: {len(pairs)} in {len(set(p[0] for p in pairs))} genomes")


if __name__ == "__main__":
    main()
