#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_annotation.py — parse multi-DB BLASTp outfmt-6 files into a structured
summary table and a homoeolog/paralog copy table (from the self-BLAST DB).

Usage:
  python summarize_annotation.py --blastdir out/blast_out --outdir out/ \
      [--self-db SelfProteome] [--top 5] [--genes GeneA GeneB]

Outputs:
  candidate_homolog_annotation_summary.csv   all top hits per gene x database
  homoeolog_copies.csv                       hits against the self proteome

Interpretation rules (apply when reading the tables; see references/pitfalls.md):
  - identity >=85% + full cover in close relatives: confident ortholog annotation
  - identity 30-60% + fragmented hits across MANY genes: repetitive/TE-derived
    protein suspected -> deprioritize cloning
  - self-BLAST hit count and identity spread = homoeolog copy inventory;
    missing copies or truncated alignments indicate CNV/annotation gaps
"""
import argparse
import csv
import sys
from pathlib import Path

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"]


def top_hits(dbfile, per_query=5):
    rows = []
    with open(dbfile) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == len(COLS):
                rows.append(dict(zip(COLS, p)))
    rows.sort(key=lambda r: (r["qseqid"], float(r["evalue"]), -float(r["bitscore"])))
    best = {}
    for r in rows:
        q = r["qseqid"].split("|")[0]
        seen = {b["sseqid"] for b in best.setdefault(q, [])}
        if r["sseqid"] not in seen:
            best[q].append(r)
    return {q: v[:per_query] for q, v in best.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blastdir", required=True, help="dir with <db>.blast.tsv files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--self-db", default="SelfProteome",
                    help="db name of the self proteome (homoeolog copies)")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--genes", nargs="*", default=None, help="restrict to these gene IDs")
    a = ap.parse_args()

    blastdir, outdir = Path(a.blastdir), Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = sorted(blastdir.glob("*.blast.tsv"))
    if not files:
        sys.exit(f"ERROR: no *.blast.tsv in {blastdir}")

    all_hits = {f.name.replace(".blast.tsv", ""): top_hits(f, a.top) for f in files}
    genes = a.genes or sorted({q for hits in all_hits.values() for q in hits})

    with open(outdir / "candidate_homolog_annotation_summary.csv", "w",
              newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["gene", "database", "rank", "hit_id", "identity_pct",
                    "aln_len", "query_cover_pct", "evalue", "bitscore"])
        for g in genes:
            for db, hits in all_hits.items():
                for i, h in enumerate(hits.get(g, []), 1):
                    w.writerow([g, db, i, h["sseqid"], h["pident"], h["length"],
                                h["qcovs"], h["evalue"], h["bitscore"]])

    if a.self_db in all_hits:
        with open(outdir / "homoeolog_copies.csv", "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(["query_gene", "copy_gene", "identity_pct", "aln_len",
                        "query_cover_pct", "evalue"])
            for g in genes:
                for h in all_hits[a.self_db].get(g, []):
                    w.writerow([g, h["sseqid"], h["pident"], h["length"],
                                h["qcovs"], h["evalue"]])

    for g in genes:
        print(f"\n### {g}")
        for db, hits in all_hits.items():
            h = hits.get(g, [])
            if h:
                t = h[0]
                print(f"  {db:20s} top={t['sseqid']:20s} id={t['pident']:>6s}% "
                      f"cov={t['qcovs']:>3s}% e={t['evalue']}")
    print("\nsummary tables written to", outdir)


if __name__ == "__main__":
    sys.exit(main())
