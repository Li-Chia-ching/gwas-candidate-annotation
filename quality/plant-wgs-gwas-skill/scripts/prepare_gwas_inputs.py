#!/usr/bin/env python3
"""Prepare sample-aligned PLINK/GEMMA phenotype, covariate, and FAM inputs."""

import argparse
import csv
import math
import shutil
import sys
from pathlib import Path


def read_fam(path):
    samples = []
    seen = set()
    with open(path, encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_no}: FAM requires at least 6 columns")
            fid, iid = fields[0], fields[1]
            key = (fid, iid)
            if iid in seen:
                raise ValueError(f"{path}:{line_no}: duplicate IID {iid}")
            seen.add(iid)
            samples.append((fid, iid, fields[:5]))
    return samples


def read_phenotypes(path, column):
    values = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"phenotype column {column} is missing from {path}")
        for line_no, row in enumerate(reader, 2):
            iid = (row.get("IID") or "").strip()
            if not iid:
                raise ValueError(f"{path}:{line_no}: IID is required")
            if iid in values:
                raise ValueError(f"{path}:{line_no}: duplicate IID {iid}")
            values[iid] = (row.get(column) or "").strip()
    return values


def read_pca(path, n_pc):
    values = {}
    with open(path, encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0].upper() == "FID":
                continue
            if len(fields) < 2 + n_pc:
                raise ValueError(f"{path}:{line_no}: fewer than {n_pc} PCs")
            fid, iid = fields[0], fields[1]
            pcs = fields[2:2 + n_pc]
            try:
                pcs = [str(float(x)) for x in pcs]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: PC values must be numeric") from exc
            values[iid] = (fid, pcs)
    return values


def read_extra_covariates(path, samples):
    if not path:
        return [], {}
    sample_iids = {iid for _, iid, _ in samples}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "IID" not in reader.fieldnames:
            raise ValueError(f"{path}: numeric covariate file requires an IID header")
        names = [x for x in reader.fieldnames if x not in ("FID", "IID")]
        values = {}
        for line_no, row in enumerate(reader, 2):
            iid = (row.get("IID") or "").strip()
            if iid not in sample_iids:
                raise ValueError(f"{path}:{line_no}: IID {iid} is not present in FAM")
            vals = []
            for name in names:
                raw = (row.get(name) or "").strip()
                try:
                    val = float(raw)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_no}: covariate {name} must be pre-encoded numeric"
                    ) from exc
                if not math.isfinite(val):
                    raise ValueError(f"{path}:{line_no}: covariate {name} is not finite")
                vals.append(raw)
            values[iid] = vals
    missing = [iid for _, iid, _ in samples if iid not in values]
    if missing:
        raise ValueError(f"{path}: covariates missing for IID {missing}")
    return names, values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fam", required=True)
    ap.add_argument("--phenotype-file", required=True)
    ap.add_argument("--phenotype-column", required=True)
    ap.add_argument("--pca", required=True)
    ap.add_argument("--n-pc", type=int, required=True)
    ap.add_argument("--extra-covariates", default="")
    ap.add_argument("--plink-pheno", required=True)
    ap.add_argument("--plink-covar", required=True)
    ap.add_argument("--gemma-fam", required=True)
    ap.add_argument("--gemma-covar", required=True)
    args = ap.parse_args()
    if args.n_pc < 0:
        raise ValueError("--n-pc must be non-negative")

    samples = read_fam(args.fam)
    phenotypes = read_phenotypes(args.phenotype_file, args.phenotype_column)
    fam_iids = {iid for _, iid, _ in samples}
    unknown = sorted(set(phenotypes) - fam_iids)
    if unknown:
        raise ValueError(f"phenotype IID values absent from FAM: {unknown}")
    pca = read_pca(args.pca, args.n_pc)
    missing_pc = [iid for _, iid, _ in samples if iid not in pca]
    if missing_pc:
        raise ValueError(f"PCA values missing for IID {missing_pc}")
    extra_names, extra = read_extra_covariates(args.extra_covariates or None, samples)

    with open(args.plink_pheno, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["FID", "IID", args.phenotype_column])
        for fid, iid, _ in samples:
            writer.writerow([fid, iid, phenotypes.get(iid, "NA") or "NA"])

    covar_header = ["FID", "IID"] + [f"PC{i}" for i in range(1, args.n_pc + 1)] + extra_names
    with open(args.plink_covar, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(covar_header)
        for fid, iid, _ in samples:
            writer.writerow([fid, iid] + pca[iid][1] + extra.get(iid, []))

    with open(args.gemma_fam, "w", encoding="utf-8") as f:
        for fid, iid, first_five in samples:
            value = phenotypes.get(iid, "")
            gemma_value = "-9" if value in ("", "NA", "NaN", "nan") else value
            f.write("\t".join([fid, iid] + first_five[2:] + [gemma_value]) + "\n")

    with open(args.gemma_covar, "w", encoding="utf-8") as f:
        for fid, iid, _ in samples:
            f.write("\t".join(["1"] + pca[iid][1] + extra.get(iid, [])) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.exit(f"ERROR: {exc}")
