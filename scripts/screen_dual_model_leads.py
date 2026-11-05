#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
screen_dual_model_leads.py — dual-model GWAS lead-SNP screening + window gene extraction.

Intersects two models' significant-signal files by exact (CHROM, POS), keeps SNPs
significant in BOTH models, clusters adjacent positions into loci (gap threshold),
takes the min-P SNP per locus as lead, and lists genes within a flanking window
from a GFF annotation.

Works with vendor-truncated signals (e.g. P<1e-5 only) as well as full summary
statistics; interpretation limits differ (see references/data_boundary.md).

Usage:
  python screen_dual_model_leads.py \
      --signals-a trait.MLM_signals.csv --signals-b trait.FarmCPU_signals.csv \
      --gff ref.gff[.gz] --trait MyTrait --outdir out/ \
      [--gap 100000] [--flank 100000] [--p-col auto]

Signals CSV must contain columns: SNP, CHROM, POS, and a P-value column
(auto-detected: the column containing 'P.value'/'P_value'/'p.value', else 8th col).
Annotation columns like snpEff_effect / snpEff_impact / snpEff_gene are carried
through when present.

Outputs (in --outdir):
  <trait>_dual_model_leads.csv      one row per independent locus (lead SNP)
  <trait>_dual_model_snps.csv       all dual-significant SNPs
  <trait>_window_genes.csv          genes within +/-flank of each lead
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def open_any(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, newline="")


def detect_pcol(header):
    for i, c in enumerate(header):
        cl = c.lower()
        if "p.value" in cl or cl in ("p", "pvalue", "p_value", "pval"):
            return i
    return 7  # vendor rMVP layout fallback (0-based 8th column)


def read_signals(path, pcol_arg):
    """Return {(chrom, pos): dict} with keys snp, chrom, pos, p, extras."""
    with open_any(path) as f:
        r = csv.DictReader(f)
        header = r.fieldnames
        pcol = pcol_arg if pcol_arg != "auto" else header[detect_pcol(header)]
        out = {}
        for row in r:
            try:
                pos = int(row["POS"])
                p = float(row[pcol])
            except (KeyError, ValueError):
                continue
            extras = {k: row.get(k, "") for k in
                      ("snpEff_effect", "snpEff_impact", "snpEff_gene", "hgvs_c", "hgvs_p",
                       "Effect", "SE") if k in row}
            out[(row["CHROM"], pos)] = dict(snp=row.get("SNP", f"{row['CHROM']}-{pos}"),
                                            chrom=row["CHROM"], pos=pos, p=p, **extras)
    return out


def load_genes(gff_path):
    genes = defaultdict(list)  # chrom -> [(start, end, gene_id)]
    with open_any(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            gid = attrs.get("ID", "")
            if gid:
                genes[p[0]].append((int(p[3]), int(p[4]), gid))
    for v in genes.values():
        v.sort()
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signals-a", required=True, help="signals CSV of model A (e.g. MLM)")
    ap.add_argument("--signals-b", required=True, help="signals CSV of model B (e.g. FarmCPU)")
    ap.add_argument("--model-a", default="A")
    ap.add_argument("--model-b", default="B")
    ap.add_argument("--gff", required=True, help="reference GFF/GFF3, optionally gzipped")
    ap.add_argument("--trait", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gap", type=int, default=100_000,
                    help="cluster gap in bp: adjacent dual SNPs closer than this merge into one locus")
    ap.add_argument("--flank", type=int, default=100_000,
                    help="window half-width for candidate gene listing")
    ap.add_argument("--p-col", default="auto", help="P-value column name (default: auto-detect)")
    a = ap.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sig_a = read_signals(a.signals_a, a.p_col)
    sig_b = read_signals(a.signals_b, a.p_col)
    dual_keys = sorted(set(sig_a) & set(sig_b))
    print(f"[{a.trait}] {a.model_a}={len(sig_a)} {a.model_b}={len(sig_b)} dual={len(dual_keys)}")

    # cluster by gap per chromosome
    by_chrom = defaultdict(list)
    for k in dual_keys:
        by_chrom[k[0]].append(k[1])
    loci = []  # (chrom, start, end, [positions])
    for chrom, poss in by_chrom.items():
        poss.sort()
        cur = [poss[0]]
        for p in poss[1:]:
            if p - cur[-1] > a.gap:
                loci.append((chrom, cur[0], cur[-1], cur))
                cur = [p]
        loci.append((chrom, cur[0], cur[-1], cur))

    # lead = min P (model A) within locus
    leads = []
    for chrom, s, e, poss in loci:
        lead_pos = min(poss, key=lambda p: min(sig_a[(chrom, p)]["p"], sig_b[(chrom, p)]["p"]))
        leads.append((chrom, lead_pos, s, e, len(poss)))
    print(f"[{a.trait}] independent loci (gap={a.gap}): {len(leads)}")

    genes = load_genes(a.gff)

    def window_genes(chrom, pos):
        return [(g, st, en) for st, en, g in genes.get(chrom, []) if en >= pos - a.flank and st <= pos + a.flank]

    with open(outdir / f"{a.trait}_dual_model_snps.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["trait", "chrom", "pos", "snp", f"P_{a.model_a}", f"P_{a.model_b}",
                    "effect_ann", "impact_ann", "gene_ann"])
        for k in dual_keys:
            ra, rb = sig_a[k], sig_b[k]
            w.writerow([a.trait, k[0], k[1], ra["snp"], ra["p"], rb["p"],
                        ra.get("snpEff_effect", "") or rb.get("snpEff_effect", ""),
                        ra.get("snpEff_impact", "") or rb.get("snpEff_impact", ""),
                        ra.get("snpEff_gene", "") or rb.get("snpEff_gene", "")])

    with open(outdir / f"{a.trait}_dual_model_leads.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["trait", "chrom", "lead_pos", "lead_snp", f"P_{a.model_a}", f"P_{a.model_b}",
                    "cluster_start", "cluster_end", "n_dual_snps", "n_window_genes",
                    "effect_ann", "impact_ann", "gene_ann"])
        for chrom, pos, s, e, n in leads:
            ra, rb = sig_a[(chrom, pos)], sig_b[(chrom, pos)]
            wg = window_genes(chrom, pos)
            w.writerow([a.trait, chrom, pos, ra["snp"], ra["p"], rb["p"], s, e, n, len(wg),
                        ra.get("snpEff_effect", "") or rb.get("snpEff_effect", ""),
                        ra.get("snpEff_impact", "") or rb.get("snpEff_impact", ""),
                        ra.get("snpEff_gene", "") or rb.get("snpEff_gene", "")])

    with open(outdir / f"{a.trait}_window_genes.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["trait", "lead_pos", "chrom", "gene_id", "gene_start", "gene_end",
                    "dist_to_lead_bp", "overlaps_lead"])
        for chrom, pos, s, e, n in leads:
            for g, st, en in window_genes(chrom, pos):
                dist = 0 if st <= pos <= en else min(abs(st - pos), abs(en - pos))
                w.writerow([a.trait, pos, chrom, g, st, en, dist, int(st <= pos <= en)])
    print(f"[{a.trait}] tables written to {outdir}")


if __name__ == "__main__":
    sys.exit(main())
