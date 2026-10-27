#!/usr/bin/env python3
import argparse
import os
import sys
import time

from Bio import Entrez


def existing_ids(path):
    ids = set()
    if not os.path.exists(path):
        return ids
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser(description="Fetch protein FASTA from NCBI by accession, resumable")
    ap.add_argument("-i", "--id-file", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--delay", type=float, default=0.4)
    args = ap.parse_args()

    if not args.email:
        sys.exit("set NCBI_EMAIL or pass --email")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    wanted = [line.strip() for line in open(args.id_file) if line.strip()]
    done = existing_ids(args.output)
    todo = [x for x in wanted if x not in done]
    print(f"total={len(wanted)} done={len(done)} todo={len(todo)}")

    with open(args.output, "a") as out:
        for i in range(0, len(todo), args.batch_size):
            batch = todo[i:i + args.batch_size]
            handle = Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text")
            text = handle.read()
            handle.close()
            if not text.endswith("\n"):
                text += "\n"
            out.write(text)
            out.flush()
            print(f"batch {i // args.batch_size + 1}: {len(batch)} records")
            time.sleep(args.delay)


if __name__ == "__main__":
    main()
