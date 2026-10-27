# MNU defense system: analysis scripts

Scripts for two parts of the paper: mining novel defense proteins in *Acinetobacter baumannii*, and tracing the Mu protein and the Mu-Nu (MNU) module across bacteria. Scripts are numbered in run order; paths are passed as arguments.

## Setup

`environment.yml` covers dRep 3.7.1, MMseqs2 18.8cc5c, BLAST+ 2.13.0, MAFFT 7.525, IQ-TREE 3.0.1 and DefensePredictor 1.1.2.

DefenseFinder 2.0.1 and PADLOC 2.0.0 run in their own environments. `06_annotate_known.sh` calls them through `conda run -n defensefinder` and `conda run -n padloc`. InterProScan 5.59-91.0 is installed separately. The download step needs the NCBI `datasets` CLI, which is not in `environment.yml`.

## Part 1: defense protein discovery

```
01  download RefSeq assemblies (chromosome level or higher)
02  dRep dereplication (fastANI, -pa 0.9 -sa 0.99 -nc 0.1)
03  DefensePredictor per genome
04  collect proteins with log-odds >= 4
05  MMseqs2 cluster at 80% identity
06  DefenseFinder + PADLOC on every genome
07  MMseqs2 search against the E. coli protein set
08  merge clusters, known systems and E. coli hits
09  InterProScan on candidate proteins
```

For example:

```bash
bash scripts/02_dereplicate_drep.sh genomes/assemblies/ncbi_dataset/data drep_out 32
python scripts/04_collect_defense_proteins.py -d rep_data -t 4.0 -o defense
python scripts/08_select_candidate_clusters.py -c defense_clusters_cluster.tsv -d rep_data -o candidate_clusters.tsv
```

`08` writes the `status`, `known_systems`, `ecoli_homolog` and `category` fields used to count known systems, novel clusters with an E. coli homolog, and novel clusters without one.

## Part 2: Mu homologs and the MNU module

```
10  PSI-BLAST of Mu against nr
11  dereplication (MMseqs2), alignment (MAFFT), tree (IQ-TREE)
12  fetch hit sequences from NCBI
13  BLASTP of Mu and Nu against a genome set
14  Mu-Nu adjacency from GFF (<= 10 kb)
15  assign GTDB R226 phylum
16  phylum distribution
```

```bash
bash scripts/10_psiblast_mu.sh MU.faa /path/to/nr 32
bash scripts/11_cluster_align_tree.sh mu_hits.faa mu 32 "LG+G4" "--fast"
python scripts/13_search_mnu.py -g refseq_genomes -m MU.faa -n NU.faa -o MNU_hits.csv
python scripts/14_check_adjacency.py -r MNU_hits.csv -g refseq_genomes -m <MU_accession> -o adjacent_MU_NU.csv
```

## Genomic neighborhoods of representative proteins

`17_extract_neighborhoods.py` takes each representative protein, locates it in its genome from the GFF, and pulls the coding genes within a nucleotide window on each side. It writes one protein FASTA per target plus a summary table. The FASTA files go straight into InterProScan.

```bash
python scripts/17_extract_neighborhoods.py \
    -t targets.tsv -g refseq_genomes -w 5000 -o neighborhoods

cat neighborhoods/*_neighborhood.faa > neighborhood_proteins.faa
bash scripts/09_interproscan.sh neighborhood_proteins.faa interproscan.tsv 16
```

`targets.tsv` has two columns, `protein_id` and `genome` (the genome directory name). Each genome directory should contain `*_genomic.gff.gz`, `*_genomic.fna.gz` and `*_protein.faa`.

## Notes

- Expected output on the A. baumannii set: 1,785 defense protein clusters, 1,065 novel candidates, 506 with an E. coli homolog (E < 1e-5).
- dRep was run with `--ignoreGenomeQuality`.
- `12_fetch_protein_fasta.py` reads `NCBI_EMAIL` and `NCBI_API_KEY` from the environment.
- `09_interproscan.sh` runs InterProScan with `-appl Pfam` by default. Pass a different analysis list as the last argument if needed.
- `11_cluster_align_tree.sh` takes a substitution model and extra IQ-TREE flags as optional arguments. The default ModelFinder search can be slow on large alignments; `--fast` shortens it.
- Large inputs (genomes, nr, GTDB R226, E. coli proteins) are not included.
