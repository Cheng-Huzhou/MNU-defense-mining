#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor


def ensure_db(genome_dir):
    faa = glob.glob(os.path.join(genome_dir, "*_protein.faa"))
    if not faa:
        return None
    db = os.path.join(genome_dir, "protein_db")
    if not os.path.exists(db + ".phr"):
        subprocess.run(["makeblastdb", "-in", faa[0], "-dbtype", "prot", "-out", db],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return db


def run_one(task):
    genome_dir, db, query_name, query_fa, evalue = task
    genome = os.path.basename(genome_dir.rstrip("/"))
    result = subprocess.run(
        ["blastp", "-query", query_fa, "-db", db, "-evalue", str(evalue),
         "-max_target_seqs", "5", "-num_threads", "1",
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore qlen slen"],
        check=True, capture_output=True, text=True)
    rows = []
    for line in result.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 8:
            continue
        rows.append([genome, query_name, p[0], p[1], p[2], p[4], p[5]])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--genome-dir", required=True)
    ap.add_argument("-m", "--mu", required=True)
    ap.add_argument("-n", "--nu", required=True)
    ap.add_argument("-e", "--evalue", type=float, default=1e-5)
    ap.add_argument("-t", "--threads", type=int, default=32)
    ap.add_argument("-o", "--output", default="MNU_hits.csv")
    args = ap.parse_args()

    genomes = sorted(d for d in glob.glob(os.path.join(args.genome_dir, "*")) if os.path.isdir(d))
    dbs = {g: ensure_db(g) for g in genomes}
    dbs = {g: db for g, db in dbs.items() if db}
    tasks = [(g, db, "MU", args.mu, args.evalue) for g, db in dbs.items()] + \
            [(g, db, "NU", args.nu, args.evalue) for g, db in dbs.items()]

    rows = []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        for res in ex.map(run_one, tasks):
            rows.extend(res)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome", "query", "qseqid", "subject", "pident", "evalue", "bitscore"])
        w.writerows(rows)

    print(f"hits: {len(rows)}")


if __name__ == "__main__":
    main()
