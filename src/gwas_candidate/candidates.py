"""GFF candidate-gene adapters: window, direct annotation, and gene-body overlap."""

import gzip
import re
from pathlib import Path

from .schemas import CandidateGeneRecord, GeneRecord, LocusRecord


def open_any(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(
        path, encoding="utf-8"
    )


def load_gff_genes(gff_path):
    genes = {}
    with open_any(gff_path) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv)
            gene_id = attrs.get("ID")
            if not gene_id:
                raise ValueError(f"{gff_path}:{line_no}: gene feature lacks ID")
            genes[gene_id] = GeneRecord(
                gene_id=gene_id,
                chrom=fields[0],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=fields[6],
            )
    return genes


def candidate_genes_for_locus(locus: LocusRecord, genes, flank):
    """Return one row per evidence source; fixed windows are not LD intervals."""
    lead = locus.lead
    locus_id = f"{locus.trait}:{locus.chrom}:{locus.start}-{locus.end}"
    records = []
    seen = set()

    def add(gene, source):
        key = (gene.gene_id, source)
        if key in seen:
            return
        seen.add(key)
        distance = 0 if gene.start <= lead.pos <= gene.end else min(
            abs(gene.start - lead.pos), abs(gene.end - lead.pos)
        )
        records.append(CandidateGeneRecord(
            trait=locus.trait,
            locus_id=locus_id,
            chrom=gene.chrom,
            lead_pos=lead.pos,
            gene_id=gene.gene_id,
            gene_start=gene.start,
            gene_end=gene.end,
            strand=gene.strand,
            distance_to_lead_bp=distance,
            overlaps_lead=gene.start <= lead.pos <= gene.end,
            evidence_source=source,
        ))

    direct = lead.annotations.get("gene_ann", "")
    for gene_id in [x for x in re.split(r"[,;|]", direct) if x]:
        if gene_id not in genes:
            raise ValueError(f"directly annotated gene {gene_id} is absent from the GFF")
        add(genes[gene_id], "direct_annotation")

    for gene in genes.values():
        if gene.chrom != locus.chrom:
            continue
        if gene.end < lead.pos - flank or gene.start > lead.pos + flank:
            continue
        if gene.start <= lead.pos <= gene.end:
            add(gene, "gene_body_overlap")
        else:
            add(gene, "fixed_window")
    source_order = {"direct_annotation": 0, "gene_body_overlap": 1, "fixed_window": 2}
    records.sort(key=lambda r: (source_order[r.evidence_source], r.distance_to_lead_bp, r.gene_id))
    return records


def candidate_genes_for_loci(loci, genes, flank):
    records = []
    for locus in loci:
        records.extend(candidate_genes_for_locus(locus, genes, flank))
    return records
