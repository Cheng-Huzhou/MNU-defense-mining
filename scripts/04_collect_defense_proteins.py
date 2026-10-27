#!/usr/bin/env python3
import argparse
import csv
import glob
import os


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--data-dir", required=True)
    ap.add_argument("-t", "--threshold", type=float, default=4.0)
    ap.add_argument("-o", "--out-prefix", default="defense")
    args = ap.parse_args()

    meta = []
    fasta = []
    for csv_path in sorted(glob.glob(os.path.join(args.data_dir, "*", "*_defense.csv"))):
        genome = os.path.basename(os.path.dirname(csv_path))
        faa = glob.glob(os.path.join(os.path.dirname(csv_path), "*_protein.faa"))
        seqs = read_fasta(faa[0]) if faa else {}

        with open(csv_path, newline="") as f:
            for row in csv.DictReader(f):
                try:
                    score = float(row["mean_log_odds"])
                except (KeyError, ValueError):
                    continue
                if score < args.threshold:
                    continue
                pid = row.get("product_accession") or row.get("non-redundant_refseq")
                if not pid or pid not in seqs:
                    continue
                meta.append([genome, pid, row.get("name", ""), row.get("locus_tag", ""),
                             row.get("start", ""), row.get("end", ""), row.get("strand", ""),
                             score])
                fasta.append((pid, seqs[pid]))

    with open(f"{args.out_prefix}_proteins.faa", "w") as f:
        for pid, seq in fasta:
            f.write(f">{pid}\n{seq}\n")

    with open(f"{args.out_prefix}_proteins_meta.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["genome", "protein_id", "name", "locus_tag", "start", "end", "strand", "mean_log_odds"])
        w.writerows(meta)

    print(f"proteins: {len(fasta)}")


if __name__ == "__main__":
    main()
