#!/usr/bin/env python3
"""Shared association-result column adapter for PLINK2, GEMMA, and standard tables."""

import csv
from dataclasses import dataclass


@dataclass(frozen=True)
class AssociationColumns:
    chrom: str
    pos: str
    p: str
    n: str | None = None
    n_missing: str | None = None
    effect: str | None = None
    se: str | None = None
    ref: str | None = None
    alt: str | None = None


def _normalized_columns(rows):
    if not rows:
        return {}
    out = {}
    for col in rows[0].keys():
        key = col.strip().lstrip("#").upper()
        out[key] = col
    return out


def _pick(cols, *names):
    for name in names:
        key = name.upper()
        if key in cols:
            return cols[key]
    return None


def read_association(path):
    """Read a delimited association file and return rows plus standardized column names."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        rows = list(csv.DictReader(f, delimiter=delimiter))
    if not rows:
        raise SystemExit("association result is empty: " + str(path))
    cols = _normalized_columns(rows)
    chrom = _pick(cols, "CHR", "CHROM", "CHROMOSOME")
    pos = _pick(cols, "BP", "POS", "POSITION", "PS")
    p = _pick(cols, "P", "PVALUE", "P_VALUE", "P_WALD")
    missing = [label for label, value in (("chromosome", chrom), ("position", pos), ("p", p)) if not value]
    if missing:
        raise SystemExit(
            "missing required association columns "
            + ", ".join(missing)
            + "; actual columns: "
            + ", ".join(rows[0].keys())
        )
    columns = AssociationColumns(
        chrom=chrom,
        pos=pos,
        p=p,
        n=_pick(cols, "N", "OBS_CT"),
        n_missing=_pick(cols, "NMISS", "N_MISS"),
        effect=_pick(cols, "BETA", "OR", "EFFECT"),
        se=_pick(cols, "SE", "STDERR", "STD_ERROR"),
        ref=_pick(cols, "REF", "ALLELE0", "A2"),
        alt=_pick(cols, "ALT", "ALLELE1", "A1"),
    )
    return rows, columns
