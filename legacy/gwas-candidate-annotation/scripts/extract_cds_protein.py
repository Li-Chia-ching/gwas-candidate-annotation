#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_cds_protein.py — extract CDS and translated proteins for target genes
from a genome FASTA + GFF, optionally building the whole proteome (longest mRNA
per gene) for homoeolog/paralog detection by self-BLAST.

Usage:
  python extract_cds_protein.py --gff ref.gff[.gz] --genome genome.fa[.gz] \
      --genes GeneA GeneB --outdir out/ [--proteome]

Outputs (in --outdir):
  candidate_cds.fasta / candidate_proteins.fasta
  extraction_report.tsv           structure sanity checks (start/stop codon,
                                  internal stops, length % 3)
  db/<genome_stem>_proteome.fasta (only with --proteome)

Sanity rules enforced and reported: ATG start, stop codon present, no internal
stop, CDS length divisible by 3. Failing any of these means wrong strand/frame
handling upstream — fix before BLAST.
"""
import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def open_any(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(seq):
    seq = seq.upper()
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--genes", nargs="+", required=True, help="gene IDs (GFF gene ID=)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--proteome", action="store_true",
                    help="also write full proteome (longest mRNA per gene) into outdir/db/")
    a = ap.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_meta, mrna_parent, cds = {}, {}, defaultdict(list)
    with open_any(a.gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            if p[2] == "gene":
                gene_meta[attrs["ID"]] = (p[0], int(p[3]), int(p[4]), p[6])
            elif p[2] in ("mRNA", "transcript"):
                if "ID" in attrs and "Parent" in attrs:
                    mrna_parent[attrs["ID"]] = attrs["Parent"]
            elif p[2] == "CDS" and "Parent" in attrs:
                cds[attrs["Parent"]].append((int(p[3]), int(p[4])))

    missing = [g for g in a.genes if g not in gene_meta]
    if missing:
        sys.exit(f"ERROR: genes not found in GFF: {missing}")

    gene_mrna = defaultdict(list)
    for m, g in mrna_parent.items():
        if g in gene_meta:
            gene_mrna[g].append(m)
    best_mrna = {g: max(ms, key=lambda m: sum(e - s + 1 for s, e in cds.get(m, [(0, 0)])))
                 for g, ms in gene_mrna.items()}

    keep = set(a.genes)
    if a.proteome:
        keep = set(gene_meta)
    chroms_needed = {gene_meta[g][0] for g in keep if g in gene_mrna}

    genome, name, buf = {}, None, []
    with open_any(a.genome) as f:
        for line in f:
            if line.startswith(">"):
                if name in chroms_needed:
                    genome[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
        if name in chroms_needed:
            genome[name] = "".join(buf)

    records, report = [], []
    for g in keep:
        if g not in best_mrna:
            continue
        m = best_mrna[g]
        chrom, gs, ge, strand = gene_meta[g]
        seq_chrom = genome.get(chrom)
        if seq_chrom is None:
            continue
        seq = "".join(seq_chrom[s - 1:e] for s, e in sorted(cds[m]))
        if strand == "-":
            seq = revcomp(seq)
        pep = translate(seq)
        records.append((g, m, chrom, gs, ge, strand, seq, pep))
        report.append(dict(gene=g, mrna=m, chrom=chrom, start=gs, end=ge, strand=strand,
                           n_cds=len(cds[m]), cds_len=len(seq), pep_len=len(pep),
                           start_codon=seq[:3], stop_codon=seq[-3:],
                           internal_stop=pep[:-1].count("*"), multi_of_3=len(seq) % 3 == 0))

    with open(outdir / "candidate_cds.fasta", "w") as f:
        for r in records:
            if r[0] in a.genes:
                f.write(f">{r[0]}|{r[1]}|{r[2]}:{r[3]}-{r[4]}({r[5]})\n{r[6]}\n")
    with open(outdir / "candidate_proteins.fasta", "w") as f:
        for r in records:
            if r[0] in a.genes:
                f.write(f">{r[0]}|{r[2]}:{r[3]}-{r[4]}({r[5]}) len={len(r[7])}aa\n{r[7]}\n")
    if a.proteome:
        (outdir / "db").mkdir(exist_ok=True)
        stem = Path(a.genome).name.replace(".fasta.gz", "").replace(".fa.gz", "").replace(".fasta", "").replace(".fa", "")
        with open(outdir / "db" / f"{stem}_proteome.fasta", "w") as f:
            for r in records:
                f.write(f">{r[0]}\n{r[7]}\n")

    keys = ["gene", "mrna", "chrom", "start", "end", "strand", "n_cds", "cds_len",
            "pep_len", "start_codon", "stop_codon", "internal_stop", "multi_of_3"]
    with open(outdir / "extraction_report.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in report:
            if r["gene"] in a.genes:
                f.write("\t".join(str(r[k]) for k in keys) + "\n")
                bad = [k for k, ok in [("start_codon", r["start_codon"] == "ATG"),
                                       ("stop", r["stop_codon"] in ("TAA", "TAG", "TGA")),
                                       ("internal_stop", r["internal_stop"] == 0),
                                       ("multi_of_3", r["multi_of_3"])] if not ok]
                status = "OK" if not bad else "CHECK:" + ",".join(bad)
                print(f"{r['gene']}: {r['pep_len']} aa, {r['n_cds']} CDS, {status}")
    print(f"done. {sum(1 for r in records if r[0] in a.genes)} target genes extracted; "
          f"{len(records)} total proteins" + (" (proteome written)" if a.proteome else ""))


if __name__ == "__main__":
    sys.exit(main())
