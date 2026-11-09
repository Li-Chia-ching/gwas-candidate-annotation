#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
screen_dual_model_leads.py — dual-model GWAS lead-SNP screening + window gene extraction.

Filters two models' association files by explicit significance thresholds, intersects
variants by chromosome, position, REF and ALT, groups adjacent positions into
distance-based clusters, selects the lead SNP from primary model A, and lists genes
within a flanking window from a GFF annotation. Distance clusters are not LD-defined
loci.

Usage:
  python screen_dual_model_leads.py \
      --signals-a trait.MLM_signals.csv --signals-b trait.FarmCPU_signals.csv \
      --gff ref.gff[.gz] --trait MyTrait --outdir out/ \
      --p-threshold-a 5e-8 --p-threshold-b 5e-8 \
      --reference-build BuildA --gff-reference-build BuildA \
      --data-mode full \
      [--gap 100000] [--flank 100000] [--p-col auto]

Signals CSV must contain SNP, CHROM, POS, REF, ALT, an effect column (Effect or BETA),
and a P-value column (auto-detected from P.value/P/pvalue unless --p-col is given).

Outputs (in --outdir):
  <trait>_dual_model_leads.csv      one row per distance cluster (lead SNP)
  <trait>_dual_model_snps.csv       all dual-significant SNPs
  <trait>_window_genes.csv          genes within +/-flank of each lead
"""
import argparse
import csv
import gzip
import math
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


def detect_named_col(header, names, label):
    lookup = {c.lower(): c for c in header}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    raise ValueError(f"missing required {label} column; available columns: {header}")


def read_signals(path, pcol_arg, p_threshold):
    """Return significant, allele-validated signals keyed by (chrom, pos, ref, alt)."""
    with open_any(path) as f:
        r = csv.DictReader(f)
        header = r.fieldnames or []
        pcol = pcol_arg if pcol_arg != "auto" else header[detect_pcol(header)]
        ref_col = detect_named_col(header, ("REF", "REF_ALLELE"), "REF")
        alt_col = detect_named_col(header, ("ALT", "ALT_ALLELE"), "ALT")
        effect_col = detect_named_col(header, ("Effect", "BETA"), "effect")
        out = {}
        coordinate_alleles = {}
        for line_no, row in enumerate(r, 2):
            try:
                chrom = row["CHROM"].strip()
                pos = int(row["POS"])
                p = float(row[pcol])
            except (KeyError, ValueError) as exc:
                raise ValueError(f"{path}:{line_no}: invalid CHROM/POS/P record") from exc
            if not (0 < p <= p_threshold):
                continue
            ref = (row.get(ref_col) or "").strip().upper()
            alt = (row.get(alt_col) or "").strip().upper()
            if ref not in ("A", "C", "G", "T") or alt not in ("A", "C", "G", "T") or ref == alt:
                raise ValueError(f"{path}:{line_no}: REF/ALT must be distinct A/C/G/T alleles")
            try:
                effect = float(row[effect_col])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_no}: effect must be numeric") from exc
            if not math.isfinite(effect) or effect == 0:
                raise ValueError(f"{path}:{line_no}: significant records require a finite non-zero effect")
            key = (chrom, pos, ref, alt)
            coord = (chrom, pos)
            previous = coordinate_alleles.get(coord)
            if previous is not None and previous != (ref, alt):
                raise ValueError(
                    f"{path}:{line_no}: multiallelic or conflicting alleles at {chrom}:{pos}; "
                    "split or normalize variants before screening"
                )
            if key in out:
                raise ValueError(f"{path}:{line_no}: duplicate variant record {chrom}:{pos} {ref}/{alt}")
            coordinate_alleles[coord] = (ref, alt)
            extras = {k: row.get(k, "") for k in
                      ("snpEff_effect", "snpEff_impact", "snpEff_gene", "hgvs_c", "hgvs_p", "SE")
                      if k in row}
            out[key] = dict(
                snp=row.get("SNP", f"{chrom}-{pos}-{ref}-{alt}"), chrom=chrom, pos=pos,
                ref=ref, alt=alt, p=p, effect=effect, **extras
            )
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
    ap.add_argument("--signals-a", required=True, help="signals CSV of primary model A")
    ap.add_argument("--signals-b", required=True, help="signals CSV of supporting model B")
    ap.add_argument("--model-a", default="A")
    ap.add_argument("--model-b", default="B")
    ap.add_argument("--gff", required=True, help="reference GFF/GFF3, optionally gzipped")
    ap.add_argument("--trait", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gap", type=int, default=100_000,
                    help="distance-cluster gap in bp; not an LD-locus definition")
    ap.add_argument("--flank", type=int, default=100_000,
                    help="window half-width for candidate gene listing")
    ap.add_argument("--p-col", default="auto", help="P-value column name (default: auto-detect)")
    ap.add_argument("--p-threshold-a", type=float, required=True,
                    help="Significance threshold applied to model A before intersection")
    ap.add_argument("--p-threshold-b", type=float, required=True,
                    help="Significance threshold applied to model B before intersection")
    ap.add_argument("--reference-build", required=True,
                    help="Reference build asserted for both association inputs")
    ap.add_argument("--gff-reference-build", required=True,
                    help="Reference build asserted for the GFF")
    ap.add_argument("--data-mode", choices=("full", "truncated"), required=True,
                    help="Explicitly declare whether summary statistics are complete or threshold-truncated")
    a = ap.parse_args()
    if not (0 < a.p_threshold_a <= 1) or not (0 < a.p_threshold_b <= 1):
        sys.exit("ERROR: P-value thresholds must be within (0, 1]")
    if a.reference_build != a.gff_reference_build:
        sys.exit(
            "ERROR: reference build mismatch: "
            f"association={a.reference_build}, GFF={a.gff_reference_build}"
        )

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sig_a = read_signals(a.signals_a, a.p_col, a.p_threshold_a)
    sig_b = read_signals(a.signals_b, a.p_col, a.p_threshold_b)

    alleles_a = {(k[0], k[1]): (k[2], k[3]) for k in sig_a}
    alleles_b = {(k[0], k[1]): (k[2], k[3]) for k in sig_b}
    for coord in sorted(set(alleles_a) & set(alleles_b)):
        if alleles_a[coord] != alleles_b[coord]:
            sys.exit(
                "ERROR: allele mismatch at "
                f"{coord[0]}:{coord[1]} ({alleles_a[coord]} vs {alleles_b[coord]})"
            )

    dual_keys = sorted(set(sig_a) & set(sig_b), key=lambda k: (k[0], k[1], k[2], k[3]))
    for key in dual_keys:
        ea, eb = sig_a[key]["effect"], sig_b[key]["effect"]
        if ea * eb < 0:
            sys.exit(
                "ERROR: opposite effect directions at "
                f"{key[0]}:{key[1]} {key[2]}/{key[3]} ({ea} vs {eb})"
            )
    print(f"[{a.trait}] {a.model_a}={len(sig_a)} {a.model_b}={len(sig_b)} dual={len(dual_keys)}")

    # Distance-based clusters are not LD-defined loci.
    by_chrom = defaultdict(list)
    for k in dual_keys:
        by_chrom[k[0]].append(k)
    clusters = []  # (chrom, start, end, [variant keys])
    for chrom, keys in by_chrom.items():
        keys.sort(key=lambda k: k[1])
        cur = [keys[0]]
        for key in keys[1:]:
            if key[1] - cur[-1][1] > a.gap:
                clusters.append((chrom, cur[0][1], cur[-1][1], cur))
                cur = [key]
            else:
                cur.append(key)
        clusters.append((chrom, cur[0][1], cur[-1][1], cur))

    # Model A is the primary analysis model; model B is cross-model support.
    leads = []
    for chrom, s, e, keys in clusters:
        lead_key = min(keys, key=lambda k: sig_a[k]["p"])
        leads.append((chrom, lead_key, s, e, len(keys)))
    print(f"[{a.trait}] distance clusters (gap={a.gap}; LD not evaluated): {len(leads)}")

    genes = load_genes(a.gff)

    def window_genes(chrom, pos):
        return [(g, st, en) for st, en, g in genes.get(chrom, [])
                if en >= pos - a.flank and st <= pos + a.flank]

    with open(outdir / f"{a.trait}_dual_model_snps.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow([
            "trait", "chrom", "pos", "snp", "ref", "alt", "reference_build", "data_mode",
            f"P_{a.model_a}", f"P_{a.model_b}", f"effect_{a.model_a}", f"effect_{a.model_b}",
            "effect_direction", "effect_ann", "impact_ann", "gene_ann"
        ])
        for k in dual_keys:
            ra, rb = sig_a[k], sig_b[k]
            w.writerow([
                a.trait, k[0], k[1], ra["snp"], k[2], k[3], a.reference_build, a.data_mode,
                ra["p"], rb["p"], ra["effect"], rb["effect"], "consistent",
                ra.get("snpEff_effect", "") or rb.get("snpEff_effect", ""),
                ra.get("snpEff_impact", "") or rb.get("snpEff_impact", ""),
                ra.get("snpEff_gene", "") or rb.get("snpEff_gene", "")
            ])

    with open(outdir / f"{a.trait}_dual_model_leads.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow([
            "trait", "chrom", "lead_pos", "lead_snp", "lead_model", "ref", "alt",
            "reference_build", "data_mode", "effect_direction", f"P_{a.model_a}", f"P_{a.model_b}",
            f"effect_{a.model_a}", f"effect_{a.model_b}",
            "cluster_start", "cluster_end", "n_dual_snps", "n_window_genes",
            "support_models", "n_support_models", "cluster_type", "ld_locus_status",
            "effect_ann", "impact_ann", "gene_ann"
        ])
        for chrom, key, s, e, n in leads:
            ra, rb = sig_a[key], sig_b[key]
            pos = key[1]
            wg = window_genes(chrom, pos)
            w.writerow([
                a.trait, chrom, pos, ra["snp"], a.model_a, key[2], key[3],
                a.reference_build, a.data_mode, "consistent", ra["p"], rb["p"], ra["effect"], rb["effect"],
                s, e, n, len(wg), f"{a.model_a};{a.model_b}", 2,
                "distance_cluster", "not_evaluated",
                ra.get("snpEff_effect", "") or rb.get("snpEff_effect", ""),
                ra.get("snpEff_impact", "") or rb.get("snpEff_impact", ""),
                ra.get("snpEff_gene", "") or rb.get("snpEff_gene", "")
            ])

    with open(outdir / f"{a.trait}_window_genes.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["trait", "lead_pos", "chrom", "gene_id", "gene_start", "gene_end",
                    "dist_to_lead_bp", "overlaps_lead"])
        for chrom, key, s, e, n in leads:
            pos = key[1]
            for g, st, en in window_genes(chrom, pos):
                dist = 0 if st <= pos <= en else min(abs(st - pos), abs(en - pos))
                w.writerow([a.trait, pos, chrom, g, st, en, dist, int(st <= pos <= en)])
    print(f"[{a.trait}] tables written to {outdir}")


if __name__ == "__main__":
    sys.exit(main())
