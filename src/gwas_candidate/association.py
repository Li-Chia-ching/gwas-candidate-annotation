"""Association-input adapters and dual-model distance clustering."""

import csv
import gzip
import math
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from .ld import pairwise_r2, read_vcf_region
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


def _chrom_sort_key(value):
    return (not value.isdigit(), int(value) if value.isdigit() else value)


def _region_dosages(vcf_path, chrom, start, end):
    """Dosage lookup by position; retries with chr-prefix toggled on naming mismatch."""
    samples, variants = read_vcf_region(vcf_path, chrom, start, end)
    if not variants:
        alt = chrom[3:] if chrom.lower().startswith("chr") else f"chr{chrom}"
        if alt != chrom:
            samples, variants = read_vcf_region(vcf_path, alt, start, end)
    return {pos: dosages for pos, dosages in variants}


def _lead_r2(dosage_by_pos, lead_a, lead_b):
    x = dosage_by_pos.get(lead_a.pos)
    y = dosage_by_pos.get(lead_b.pos)
    if x is None or y is None:
        return None
    value = pairwise_r2(x, y)
    return None if math.isnan(value) else value


def _merged_record(entries, indices, evaluated):
    members = [entries[i] for i in indices]
    if len(members) == 1:
        locus = members[0]
        status = "ld_evaluated_not_merged" if evaluated[indices[0]] else "not_evaluated"
        return replace(locus, ld_locus_status=status)
    variants = sorted((v for locus in members for v in locus.variants), key=lambda v: v.pos)
    lead = min(variants, key=lambda v: v.p)
    support = tuple(dict.fromkeys(m for locus in members for m in locus.support_models))
    return LocusRecord(
        trait=members[0].trait,
        chrom=members[0].chrom,
        start=min(locus.start for locus in members),
        end=max(locus.end for locus in members),
        lead=lead,
        lead_model=lead.model,
        support_models=support,
        variants=tuple(variants),
        reference_build=members[0].reference_build,
        data_mode=members[0].data_mode,
        cluster_type="ld_merged_cluster",
        ld_locus_status="merged_by_ld",
        merged_by_ld=True,
        merged_cluster_count=len(members),
    )


def merge_clusters_by_ld(loci, vcf_path, r2_threshold=0.5):
    """Merge adjacent distance clusters whose lead SNPs exceed an LD r2 threshold.

    Single-linkage over adjacent clusters on the same chromosome: when the
    lead-lead dosage r2 is above ``r2_threshold``, both clusters join one LD
    locus and the merged locus takes the minimum-P variant as its nominal lead.
    Evaluation is iterative, so a merged locus is re-tested against its next
    neighbour with the recomputed lead. The lead SNP is a nominal lead only;
    no conditional analysis of multiple independent signals is performed.

    Loci whose lead genotypes are absent from the VCF (or all loci when
    ``vcf_path`` is None) are returned unchanged with ``ld_locus_status``
    ``not_evaluated``.
    """
    if not (0 < r2_threshold <= 1):
        raise ValueError("r2_threshold must be within (0, 1]")
    entries = tuple(loci)
    if len(entries) < 2:
        return entries
    if vcf_path is None:
        return entries
    by_chrom = defaultdict(list)
    for locus in entries:
        by_chrom[locus.chrom].append(locus)
    result = []
    for chrom in sorted(by_chrom, key=_chrom_sort_key):
        group = sorted(by_chrom[chrom], key=lambda l: (l.start, l.end, l.lead.pos))
        dosage_by_pos = _region_dosages(
            vcf_path, chrom, min(l.start for l in group), max(l.end for l in group)
        )
        blocks = [[i] for i in range(len(group))]
        evaluated = [False] * len(group)
        changed = True
        while changed and len(blocks) > 1:
            changed = False
            merged_blocks = []
            i = 0
            while i < len(blocks):
                if i + 1 < len(blocks):
                    left, right = blocks[i], blocks[i + 1]
                    lead_left = min((group[j].lead for j in left), key=lambda v: v.p)
                    lead_right = min((group[j].lead for j in right), key=lambda v: v.p)
                    r2 = _lead_r2(dosage_by_pos, lead_left, lead_right)
                    if r2 is not None:
                        for j in left + right:
                            evaluated[j] = True
                    if r2 is not None and r2 > r2_threshold:
                        merged_blocks.append(left + right)
                        changed = True
                        i += 2
                        continue
                merged_blocks.append(blocks[i])
                i += 1
            blocks = merged_blocks
        result.extend(_merged_record(group, block, evaluated) for block in blocks)
    return tuple(result)


def rank_loci(loci):
    """Rank loci by nominal lead P (ascending); rank 1 is the priority peak.

    Returns loci sorted by rank with ``locus_rank`` and ``priority_peak``
    assigned. Ties break on chromosome, start, then end for determinism.
    """
    ordered = sorted(
        loci, key=lambda l: (l.lead.p, _chrom_sort_key(l.chrom), l.start, l.end)
    )
    return tuple(
        replace(locus, locus_rank=rank, priority_peak=(rank == 1))
        for rank, locus in enumerate(ordered, 1)
    )
