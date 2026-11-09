"""Pairwise LD (r2) from biallelic VCF genotype dosages."""

import gzip
import math
from pathlib import Path

import numpy as np


def open_any(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(
        path, encoding="utf-8"
    )


def read_vcf_region(path, chrom, start, end):
    variants = []
    with open_any(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] != chrom:
                continue
            pos = int(fields[1])
            if pos < start or pos > end:
                continue
            if "," in fields[4]:
                raise ValueError(f"multiallelic VCF record is not supported for LD: {chrom}:{pos}")
            fmt = fields[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError as exc:
                raise ValueError(f"GT is missing at {chrom}:{pos}") from exc
            dosages = []
            for sample in fields[9:]:
                gt = sample.split(":")[gt_index]
                if "." in gt:
                    dosages.append(None)
                    continue
                sep = "/" if "/" in gt else "|"
                try:
                    dosages.append(sum(int(x) for x in gt.split(sep)))
                except ValueError as exc:
                    raise ValueError(f"invalid GT at {chrom}:{pos}: {gt}") from exc
            variants.append((pos, dosages))
    return samples, variants


def pairwise_r2(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 2:
        return math.nan
    a = np.array([p[0] for p in pairs], dtype=float)
    b = np.array([p[1] for p in pairs], dtype=float)
    if np.var(a) == 0 or np.var(b) == 0:
        return math.nan
    return float(np.corrcoef(a, b)[0, 1] ** 2)


def ld_matrix_for_region(vcf_path, chrom, start, end):
    samples, variants = read_vcf_region(vcf_path, chrom, start, end)
    labels = [str(pos) for pos, _ in variants]
    matrix = np.full((len(variants), len(variants)), np.nan)
    for i, (_, x) in enumerate(variants):
        for j, (_, y) in enumerate(variants):
            matrix[i, j] = 1.0 if i == j else pairwise_r2(x, y)
    return labels, matrix, len(samples)


def write_ld_matrix(path, labels, matrix):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("variant\t" + "\t".join(labels) + "\n")
        for label, row in zip(labels, matrix):
            handle.write(label + "\t" + "\t".join("NA" if math.isnan(x) else f"{x:.6f}" for x in row) + "\n")
