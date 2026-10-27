#!/usr/bin/env python3
import argparse
import csv
import glob
import gzip
import os


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def parse_gff_cds(path):
    feats = []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            pid = attrs.get("protein_id")
            if not pid:
                continue
            feats.append({
                "contig": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "protein_id": pid,
                "product": attrs.get("product", ""),
            })
    return feats


def first_existing(genome_dir, patterns, exclude=None):
    for pattern in patterns:
        for path in sorted(glob.glob(os.path.join(genome_dir, pattern))):
            if exclude and exclude in os.path.basename(path):
                continue
            return path
    return None


def revcomp(seq):
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def main():
    ap = argparse.ArgumentParser(description="Extract genomic neighborhoods around target proteins")
    ap.add_argument("-t", "--targets", required=True, help="TSV with columns protein_id and genome")
    ap.add_argument("-g", "--genome-dir", required=True)
    ap.add_argument("-w", "--window", type=int, default=5000)
    ap.add_argument("-o", "--outdir", default="neighborhoods")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    summary = open(os.path.join(args.outdir, "neighborhood_summary.tsv"), "w", newline="")
    sw = csv.writer(summary, delimiter="\t")
    sw.writerow(["target", "genome", "contig", "target_start", "target_end", "strand",
                 "window_start", "window_end", "n_genes"])

    genes_out = open(os.path.join(args.outdir, "neighborhood_genes.tsv"), "w", newline="")
    gw = csv.writer(genes_out, delimiter="\t")
    gw.writerow(["target", "genome", "contig", "protein_id", "start", "end", "strand", "product"])

    with open(args.targets, newline="") as f:
        targets = list(csv.DictReader(f, delimiter="\t"))

    for row in targets:
        protein_id = row["protein_id"].strip()
        genome = row["genome"].strip()
        genome_dir = os.path.join(args.genome_dir, genome)

        gff = first_existing(genome_dir, ["*_genomic.gff.gz", "*_genomic.gff", "*.gff.gz", "*.gff"])
        faa = first_existing(genome_dir, ["*_protein.faa", "*.faa"])
        fna = first_existing(genome_dir, ["*_genomic.fna", "*_genomic.fna.gz", "*.fna", "*.fna.gz"],
                             exclude="cds_from_genomic")
        if not gff or not faa:
            print(f"skip {protein_id}: missing gff/faa in {genome_dir}")
            continue

        feats = parse_gff_cds(gff)
        target = next((x for x in feats if x["protein_id"] == protein_id), None)
        if target is None:
            print(f"skip {protein_id}: not found in {genome}")
            continue

        contig = target["contig"]
        win_start = max(0, target["start"] - args.window)
        win_end = target["end"] + args.window

        in_window = [x for x in feats
                     if x["contig"] == contig and x["start"] <= win_end and x["end"] >= win_start]
        in_window.sort(key=lambda x: x["start"])

        seqs = read_fasta(faa)
        out_faa = os.path.join(args.outdir, f"{protein_id}_neighborhood.faa")
        with open(out_faa, "w") as out:
            for gene in in_window:
                seq = seqs.get(gene["protein_id"])
                if not seq:
                    continue
                product = gene["product"].replace(" ", "_").replace("|", "_") or "hypothetical"
                out.write(f">{gene['protein_id']}|{contig}|{gene['start']}-{gene['end']}|{gene['strand']}|{product}\n")
                out.write(f"{seq}\n")
                gw.writerow([protein_id, genome, contig, gene["protein_id"],
                             gene["start"], gene["end"], gene["strand"], gene["product"]])

        if fna:
            genome_seqs = read_fasta(fna)
            contig_seq = genome_seqs.get(contig, "")
            if contig_seq:
                window_seq = contig_seq[win_start - 1:win_end]
                if target["strand"] == "-":
                    window_seq = revcomp(window_seq)
                with open(os.path.join(args.outdir, f"{protein_id}_neighborhood.fna"), "w") as out:
                    out.write(f">{contig}:{win_start}-{win_end}({target['strand']})\n{window_seq}\n")

        sw.writerow([protein_id, genome, contig, target["start"], target["end"], target["strand"],
                     win_start, win_end, len(in_window)])

    summary.close()
    genes_out.close()


if __name__ == "__main__":
    main()
