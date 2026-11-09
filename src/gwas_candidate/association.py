"""Association-input adapters and dual-model distance clustering."""

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from .schemas import LocusRecord, VariantRecord


def open_any(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix == ".gz" else open(
        path, newline="", encoding="utf-8-sig"
    )


def detect_pcol(header):
    for i, col in enumerate(header):
        lower = col.lower()
        if "p.value" in lower or lower in ("p", "pvalue", "p_value", "pval", "p_wald"):
            return i
    raise ValueError(f"missing P-value column; available columns: {header}")


def _named_col(header, names, label):
    lookup = {col.lower(): col for col in header}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    raise ValueError(f"missing required {label} column; available columns: {header}")


def read_model_signals(path, model, trait, p_threshold, data_mode, p_col="auto"):
    """Read one model's association CSV/TSV after explicit significance filtering."""
    if data_mode not in ("full", "truncated"):
        raise ValueError("data_mode must be 'full' or 'truncated'")
    if not (0 < p_threshold <= 1):
        raise ValueError("p_threshold must be within (0, 1]")
    with open_any(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = reader.fieldnames or []
        pcol = p_col if p_col != "auto" else header[detect_pcol(header)]
        chrom_col = _named_col(header, ("CHROM", "CHR", "#CHROM"), "chromosome")
        pos_col = _named_col(header, ("POS", "BP", "PS"), "position")
        ref_col = _named_col(header, ("REF", "REF_ALLELE", "ALLELE0", "A2"), "REF")
        alt_col = _named_col(header, ("ALT", "ALT_ALLELE", "ALLELE1", "A1"), "ALT")
        effect_col = _named_col(header, ("Effect", "BETA"), "effect")
        snp_col = _named_col(header, ("SNP", "ID", "RS"), "SNP")
        records = {}
        coordinate_alleles = {}
        for line_no, row in enumerate(reader, 2):
            try:
                chrom = row[chrom_col].strip()
                pos = int(float(row[pos_col]))
                p = float(row[pcol])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_no}: invalid chromosome/position/P record") from exc
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
                raise ValueError(f"{path}:{line_no}: significant records require finite non-zero effects")
            key = (chrom, pos, ref, alt)
            coord = (chrom, pos)
            previous = coordinate_alleles.get(coord)
            if previous is not None and previous != (ref, alt):
                raise ValueError(
                    f"{path}:{line_no}: multiallelic or conflicting alleles at {chrom}:{pos}"
                )
            if key in records:
                raise ValueError(f"{path}:{line_no}: duplicate variant record {chrom}:{pos} {ref}/{alt}")
            coordinate_alleles[coord] = (ref, alt)
            annotations = {
                name: row.get(source, "")
                for name, source in (
                    ("effect_ann", "snpEff_effect"),
                    ("impact_ann", "snpEff_impact"),
                    ("gene_ann", "snpEff_gene"),
                    ("hgvs_c", "hgvs_c"),
                    ("hgvs_p", "hgvs_p"),
                    ("se", "SE"),
                )
                if source in row
            }
            records[key] = VariantRecord(
                trait=trait, model=model, chrom=chrom, pos=pos,
                snp=row.get(snp_col) or f"{chrom}:{pos}:{ref}:{alt}",
                ref=ref, alt=alt, p=p, effect=effect, data_mode=data_mode,
                annotations=annotations,
            )
    return records


def build_dual_model_clusters(
    signals_a,
    signals_b,
    model_a,
    model_b,
    trait,
    gap,
    reference_build,
    data_mode,
):
    """Build distance clusters; no LD or independence claim is made."""
    if data_mode not in ("full", "truncated"):
        raise ValueError("data_mode must be 'full' or 'truncated'")
    alleles_a = {(k[0], k[1]): (k[2], k[3]) for k in signals_a}
    alleles_b = {(k[0], k[1]): (k[2], k[3]) for k in signals_b}
    for coord in sorted(set(alleles_a) & set(alleles_b)):
        if alleles_a[coord] != alleles_b[coord]:
            raise ValueError(f"allele mismatch at {coord[0]}:{coord[1]}")
    shared_keys = sorted(set(signals_a) & set(signals_b), key=lambda k: (k[0], k[1], k[2], k[3]))
    for key in shared_keys:
        if signals_a[key].effect * signals_b[key].effect < 0:
            raise ValueError(f"opposite effect directions at {key[0]}:{key[1]} {key[2]}/{key[3]}")
    by_chrom = defaultdict(list)
    for key in shared_keys:
        by_chrom[key[0]].append(key)
    loci = []
    for chrom, keys in by_chrom.items():
        keys.sort(key=lambda k: k[1])
        current = [keys[0]]
        for key in keys[1:]:
            if key[1] - current[-1][1] > gap:
                loci.append(_cluster_record(current, signals_a, signals_b, model_a, model_b, trait, reference_build, data_mode))
                current = [key]
            else:
                current.append(key)
        loci.append(_cluster_record(current, signals_a, signals_b, model_a, model_b, trait, reference_build, data_mode))
    return loci


def _cluster_record(keys, signals_a, signals_b, model_a, model_b, trait, reference_build, data_mode):
    lead_key = min(keys, key=lambda k: signals_a[k].p)
    variants = tuple(signals_a[k] for k in keys)
    return LocusRecord(
        trait=trait,
        chrom=keys[0][0],
        start=min(k[1] for k in keys),
        end=max(k[1] for k in keys),
        lead=signals_a[lead_key],
        lead_model=model_a,
        support_models=(model_a, model_b),
        variants=variants,
        reference_build=reference_build,
        data_mode=data_mode,
    )
