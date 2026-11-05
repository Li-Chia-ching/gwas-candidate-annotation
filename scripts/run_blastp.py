#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_blastp.py — multi-database BLASTp for candidate proteins (local NCBI BLAST+).

Databases are declared in a two-column TSV config:  name<TAB>/path/to/proteins.fasta
Recommended tiering (see references/workflow.md):
  1) same-species reference proteome (e.g. your crop's reference)
  2) close-relative species proteome
  3) model-plant proteome (e.g. Arabidopsis TAIR10)
  4) self proteome from extract_cds_protein.py --proteome  (homoeolog/paralog copies)

Usage:
  python run_blastp.py --query candidate_proteins.fasta --dbs dbs.tsv --outdir out/ \
      [--blastp /path/blastp] [--makeblastdb /path/makeblastdb] \
      [--evalue 1e-5] [--max-target 20] [--threads 8]

Output per DB: blast_out/<name>.blast.tsv
  outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send
            evalue bitscore qcovs
"""
import argparse
import subprocess
import sys
from pathlib import Path

OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
          "sstart send evalue bitscore qcovs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--dbs", required=True, help="TSV: db_name<TAB>protein_fasta_path")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--blastp", default="blastp")
    ap.add_argument("--makeblastdb", default="makeblastdb")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--max-target", default="20")
    ap.add_argument("--threads", default="8")
    a = ap.parse_args()

    outdir = Path(a.outdir)
    (outdir / "blast_out").mkdir(parents=True, exist_ok=True)

    dbs = []
    with open(a.dbs) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            name, path = line.split("\t")[:2]
            dbs.append((name, Path(path)))

    for name, fasta in dbs:
        if not fasta.exists():
            sys.exit(f"ERROR: DB fasta missing for {name}: {fasta}")
        dbdir = outdir / "db" / name
        dbdir.mkdir(parents=True, exist_ok=True)
        dbpfx = dbdir / name
        if not (dbdir / f"{name}.psq").exists():
            subprocess.run([a.makeblastdb, "-in", str(fasta), "-dbtype", "prot",
                            "-out", str(dbpfx), "-parse_seqids"], check=True)
            print(f"[makedb] {name}")
        out = outdir / "blast_out" / f"{name}.blast.tsv"
        subprocess.run([a.blastp, "-query", a.query, "-db", str(dbpfx),
                        "-out", str(out), "-outfmt", OUTFMT,
                        "-evalue", a.evalue, "-max_target_seqs", a.max_target,
                        "-num_threads", a.threads], check=True)
        print(f"[blastp] {name} -> {out.name}")
    print("ALL BLAST DONE")


if __name__ == "__main__":
    sys.exit(main())
