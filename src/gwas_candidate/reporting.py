"""Candidate-gene report and figure generation.

Figures use English labels, Arial, Okabe-Ito colors, and 600-dpi PDF export.
"""

import csv
import html
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .ld import ld_matrix_for_region, write_ld_matrix
from .schemas import CandidateGeneRecord, HomologyHit, LocusRecord, ProteinRecord


OKABE_BLUE = "#0072B2"
OKABE_VERMILION = "#D55E00"
OKABE_ORANGE = "#E69F00"
OKABE_GREEN = "#009E73"
OKABE_SKY = "#56B4E9"


def _style():
    plt.rcParams.update({
        "font.family": "Arial",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def _save(fig, path):
    fig.savefig(path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def read_association_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    if not rows:
        raise ValueError(f"association table is empty: {path}")
    lookup = {c.strip().lstrip("#").upper(): c for c in rows[0]}

    def pick(*names):
        for name in names:
            if name in lookup:
                return lookup[name]
        raise ValueError(f"missing association column among {names}")

    chrom_col = pick("CHR", "CHROM", "CHROMOSOME")
    pos_col = pick("BP", "POS", "POSITION", "PS")
    p_col = pick("P", "P.VALUE", "PVALUE", "P_VALUE", "P_WALD")
    parsed = []
    for row in rows:
        try:
            p = float(row[p_col])
            pos = int(float(row[pos_col]))
        except (TypeError, ValueError):
            continue
        if 0 < p <= 1:
            parsed.append({"chrom": str(row[chrom_col]).replace("chr", ""), "pos": pos, "p": p})
    return parsed


def _chrom_key(value):
    return (not value.isdigit(), int(value) if value.isdigit() else value)


def plot_manhattan(rows, path, threshold):
    _style()
    chroms = sorted({r["chrom"] for r in rows}, key=_chrom_key)
    fig, ax = plt.subplots(figsize=(8, 3.2))
    offset = 0
    ticks = []
    for i, chrom in enumerate(chroms):
        subset = [r for r in rows if r["chrom"] == chrom]
        x = np.array([r["pos"] for r in subset]) - min(r["pos"] for r in subset) + offset
        y = -np.log10([r["p"] for r in subset])
        ax.scatter(x, y, s=7, color=OKABE_BLUE if i % 2 == 0 else OKABE_SKY, linewidths=0)
        ticks.append(offset + (x.max() - x.min()) / 2)
        offset = x.max() + 1
    ax.axhline(-math.log10(threshold), color=OKABE_VERMILION, linestyle="--", linewidth=1.2)
    ax.set_xticks(ticks, chroms, fontsize=7)
    ax.set_xlabel("Genomic position (chromosome)")
    ax.set_ylabel("-log10(P)")
    _save(fig, path)


def plot_qq(rows, path):
    _style()
    p = np.sort(np.array([r["p"] for r in rows]))
    observed = -np.log10(p)[::-1]
    expected = (-np.log10(np.arange(1, len(p) + 1) / (len(p) + 1)))[::-1]
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    limit = max(expected.max(), observed.max()) * 1.05
    ax.plot([0, limit], [0, limit], color=OKABE_VERMILION, linestyle="--", linewidth=1.2)
    ax.plot(expected, observed, color=OKABE_BLUE, linewidth=1.5)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    _save(fig, path)


def plot_regional(rows, locus, candidates, path, flank):
    _style()
    chrom = locus.chrom.replace("chr", "")
    left = max(1, locus.start - flank)
    right = locus.end + flank
    subset = [r for r in rows if r["chrom"] == chrom and left <= r["pos"] <= right]
    fig, (ax, ax_gene) = plt.subplots(
        2, 1, figsize=(7, 4.2), sharex=True, height_ratios=[3, 1],
        gridspec_kw={"hspace": 0.08},
    )
    if subset:
        ax.scatter(
            [r["pos"] for r in subset], [-math.log10(r["p"]) for r in subset],
            s=14, color=OKABE_BLUE, linewidths=0,
        )
    ax.axvline(locus.lead.pos, color=OKABE_VERMILION, linestyle="--", linewidth=1.1)
    ax.set_ylabel("-log10(P)")
    ax.set_xlim(left, right)
    gene_rows = [c for c in candidates if c.locus_id == _locus_id(locus)]
    unique_genes = sorted({c.gene_id: c for c in gene_rows}.values(), key=lambda c: c.gene_start)
    for i, candidate in enumerate(unique_genes):
        color = OKABE_ORANGE if candidate.overlaps_lead else OKABE_BLUE
        ax_gene.hlines(i, candidate.gene_start, candidate.gene_end, color=color, linewidth=5)
    ax_gene.set_yticks(range(len(unique_genes)), [c.gene_id for c in unique_genes])
    ax_gene.set_xlabel(f"Chromosome {locus.chrom} position (bp)")
    ax_gene.set_ylabel("Genes")
    _save(fig, path)


def plot_ld_heatmap(labels, matrix, path):
    _style()
    cmap = LinearSegmentedColormap.from_list("okabe_ld", ["#FFFFFF", OKABE_BLUE])
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    masked = np.ma.masked_invalid(matrix)
    image = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.set_xlabel("Variant position")
    ax.set_ylabel("Variant position")
    fig.colorbar(image, ax=ax, label="r²")
    _save(fig, path)


def plot_evidence_matrix(candidates, proteins, hits, assessments, path):
    _style()
    genes = sorted({c.gene_id for c in candidates})
    columns = [
        "Direct annotation", "Gene body", "Window", "Protein QC",
        "Homology candidate", "TE/repeat risk",
    ]
    candidate_sources = defaultdict(set)
    for c in candidates:
        candidate_sources[c.gene_id].add(c.evidence_source)
    protein_ok = {p.gene_id: p.qc_status == "OK" for p in proteins}
    hit_genes = {h.query_gene for h in hits}
    matrix = []
    for gene in genes:
        assessment = assessments.get(gene)
        matrix.append([
            "direct_annotation" in candidate_sources[gene],
            "gene_body_overlap" in candidate_sources[gene],
            "fixed_window" in candidate_sources[gene],
            protein_ok.get(gene, False),
            gene in hit_genes,
            assessment is not None and assessment.te_status == "detected",
        ])
    fig_height = max(2.4, 0.35 * len(genes) + 1.4)
    fig, ax = plt.subplots(figsize=(6.2, fig_height))
    ax.imshow(np.array(matrix, dtype=float), cmap=LinearSegmentedColormap.from_list(
        "evidence", ["#FFFFFF", OKABE_GREEN]
    ), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(columns)), columns, rotation=35, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    ax.set_xlabel("Evidence class")
    ax.set_ylabel("Candidate gene")
    _save(fig, path)


def _variant_protein_position(locus, protein):
    text = locus.lead.annotations.get("hgvs_p", "")
    match = re.search(r"p\.[A-Za-z*]+(\d+)", text)
    if match:
        return int(match.group(1))
    span = max(1, locus.end - locus.start + 1)
    ratio = min(max((locus.lead.pos - protein.start) / max(1, protein.end - protein.start), 0), 1)
    return max(1, round(ratio * max(1, len(protein.protein_sequence.rstrip("*")))))


def plot_domain_maps(loci, proteins, domains, figures_dir):
    _style()
    paths = {}
    proteins_by_gene = {p.gene_id: p for p in proteins}
    lead_by_gene = {}
    for locus in loci:
        for candidate_gene in getattr(locus, "candidate_gene_ids", ()):  # optional prelink
            lead_by_gene.setdefault(candidate_gene, locus)
    for gene_id, protein in proteins_by_gene.items():
        locus = lead_by_gene.get(gene_id, loci[0])
        length = max(1, len(protein.protein_sequence.rstrip("*")))
        variant_pos = min(_variant_protein_position(locus, protein), length)
        fig, ax = plt.subplots(figsize=(6.5, 1.8))
        ax.hlines(0, 1, length, color=OKABE_BLUE, linewidth=5)
        for domain in domains.get(gene_id, []):
            start = max(1, int(domain["start"]))
            end = min(length, int(domain["end"]))
            ax.hlines(0, start, end, color=OKABE_ORANGE, linewidth=9)
            ax.text((start + end) / 2, 0.12, str(domain.get("name", "domain")), ha="center", fontsize=8)
        ax.scatter([variant_pos], [0], s=55, color=OKABE_VERMILION, marker="v", zorder=3)
        ax.text(variant_pos, -0.16, "variant", ha="center", fontsize=8)
        ax.set_xlim(0, length + 1)
        ax.set_ylim(-0.35, 0.35)
        ax.set_yticks([])
        ax.set_xlabel("Protein position (aa)")
        path = figures_dir / f"domain_map_{gene_id}.pdf"
        _save(fig, path)
        paths[gene_id] = path
    return paths


def _locus_id(locus):
    return f"{locus.trait}:{locus.chrom}:{locus.start}-{locus.end}"


def _write_tables(outdir, loci, candidates, hits, assessments):
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    with (tables / "loci.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "locus_id", "locus_rank", "priority_peak", "trait", "chrom", "start", "end",
            "lead_snp", "lead_pos", "lead_p", "lead_model", "support_models",
            "cluster_type", "ld_locus_status", "merged_by_ld", "merged_cluster_count",
            "data_mode", "reference_build",
        ])
        for locus in loci:
            writer.writerow([
                _locus_id(locus), locus.locus_rank if locus.locus_rank is not None else "",
                int(locus.priority_peak), locus.trait, locus.chrom, locus.start, locus.end,
                locus.lead.snp, locus.lead.pos, f"{locus.lead.p:.6g}", locus.lead_model,
                ";".join(locus.support_models), locus.cluster_type,
                locus.ld_locus_status, int(locus.merged_by_ld), locus.merged_cluster_count,
                locus.data_mode, locus.reference_build,
            ])
    with (tables / "candidate_genes.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "locus_id", "gene_id", "chrom", "gene_start", "gene_end", "strand",
            "distance_to_lead_bp", "overlaps_lead", "evidence_source",
        ])
        for c in candidates:
            writer.writerow([
                c.locus_id, c.gene_id, c.chrom, c.gene_start, c.gene_end, c.strand,
                c.distance_to_lead_bp, int(c.overlaps_lead), c.evidence_source,
            ])
    with (tables / "homology_candidates.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_gene", "database", "hit_id", "identity_pct", "alignment_length",
            "query_cover_pct", "subject_cover_pct", "evalue", "bitscore", "evidence_label",
        ])
        for h in hits:
            writer.writerow([
                h.query_gene, h.database, h.hit_id, h.identity_pct, h.alignment_length,
                h.query_cover_pct, h.subject_cover_pct, h.evalue, h.bitscore, h.evidence_label,
            ])
    with (tables / "candidate_assessments.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene_id", "priority", "evidence_level", "te_status", "cloning_gate",
            "hormone_rescue_gate", "recommendation", "downgrade_reasons",
        ])
        for gene_id in sorted(assessments):
            a = assessments[gene_id]
            writer.writerow([
                a.gene_id, a.priority, a.evidence_level, a.te_status, a.cloning_gate,
                a.hormone_rescue_gate, a.recommendation, ";".join(a.downgrade_reasons),
            ])
    return tables


def generate_report(
    trait,
    association_path,
    loci,
    candidates,
    proteins,
    homology_hits,
    assessments,
    domains,
    outdir,
    threshold,
    data_mode,
    reference_build,
    vcf_path=None,
    region_flank=100000,
):
    if data_mode not in ("full", "truncated"):
        raise ValueError("data_mode must be full or truncated")
    outdir = Path(outdir)
    figures = outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    rows = read_association_rows(association_path)
    manhattan = figures / f"{trait}_manhattan.pdf"
    plot_manhattan(rows, manhattan, threshold)
    qq = None
    if data_mode == "full":
        qq = figures / f"{trait}_qq.pdf"
        plot_qq(rows, qq)

    regional_paths = []
    ld_paths = []
    ld_tables = []
    for index, locus in enumerate(loci, 1):
        regional = figures / f"locus_{index:02d}_regional_association.pdf"
        plot_regional(rows, locus, candidates, regional, region_flank)
        regional_paths.append(regional)
        if vcf_path:
            labels, matrix, sample_count = ld_matrix_for_region(
                vcf_path, locus.chrom, max(1, locus.start - region_flank), locus.end + region_flank
            )
            ld_path = figures / f"locus_{index:02d}_ld_heatmap.pdf"
            plot_ld_heatmap(labels, matrix, ld_path)
            ld_paths.append(ld_path)
            table = outdir / "tables" / f"locus_{index:02d}_ld_matrix.tsv"
            table.parent.mkdir(parents=True, exist_ok=True)
            write_ld_matrix(table, labels, matrix)
            ld_tables.append(table)

    evidence_pdf = figures / "candidate_evidence_matrix.pdf"
    plot_evidence_matrix(candidates, proteins, homology_hits, assessments, evidence_pdf)
    domain_paths = plot_domain_maps(loci, proteins, domains, figures)
    tables = _write_tables(outdir, loci, candidates, homology_hits, assessments)

    assessment_rows = []
    for gene_id in sorted(assessments):
        a = assessments[gene_id]
        assessment_rows.append([
            gene_id, a.priority, a.evidence_level, a.te_status, a.cloning_gate,
            a.hormone_rescue_gate, a.recommendation, ";".join(a.downgrade_reasons) or "none",
        ])

    locus_rows = []
    for locus in loci:
        locus_rows.append([
            locus.locus_rank if locus.locus_rank is not None else "-",
            _locus_id(locus),
            f"{locus.lead.snp} ({locus.lead.pos})",
            f"{locus.lead.p:.3g}",
            locus.cluster_type,
            locus.ld_locus_status,
            "yes" if locus.merged_by_ld else "no",
        ])

    md = outdir / "report.md"
    html_path = outdir / "report.html"
    figure_lines = [f"- figures/{manhattan.name}"]
    if qq:
        figure_lines.append(f"- figures/{qq.name}")
    else:
        figure_lines.append("- Q-Q: not evaluated for truncated input")
    figure_lines += [f"- figures/{p.name}" for p in regional_paths]
    figure_lines += [f"- figures/{p.name}" for p in ld_paths]
    figure_lines.append(f"- figures/{evidence_pdf.name}")
    figure_lines += [f"- figures/{p.name}" for p in domain_paths.values()]

    md_text = "\n".join([
        f"# Candidate-gene report: {trait}",
        "",
        "V3 reports candidate Priority A/B/C/D and Evidence L0-L4 only; it does not claim causal-gene confirmation.",
        "",
        "Lead SNPs are nominal leads (minimum P within each locus); conditional analysis for multiple independent signals within one locus is not implemented.",
        "",
        f"- Data mode: `{data_mode}`",
        f"- Reference genome: `{reference_build}`",
        f"- Threshold: `{threshold:g}`",
        f"- Loci: {len(loci)}; candidate evidence rows: {len(candidates)}",
        "",
        "## Loci (ranked by nominal lead P)",
        "",
        "| Rank | Locus | Lead SNP (pos) | Lead P | Cluster type | LD locus status | Merged by LD |",
        "| --- | --- | --- | --- | --- | --- | --- |",
        *[
            "| " + " | ".join(str(x).replace("|", "/") for x in row) + " |"
            for row in locus_rows
        ],
        "",
        "## Figures",
        *figure_lines,
        "",
        "## Candidate assessments",
        "",
        "| Gene | Priority | Evidence | TE | Cloning | Hormone rescue | Recommendation | Downgrade reasons |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
        *[
            "| " + " | ".join(str(x).replace("|", "/") for x in row) + " |"
            for row in assessment_rows
        ],
        "",
        "## Standard tables",
        "- `tables/loci.tsv`",
        "- `tables/candidate_genes.tsv`",
        "- `tables/homology_candidates.tsv`",
        "- `tables/candidate_assessments.tsv`",
        "",
    ])
    md.write_text(md_text, encoding="utf-8")

    cards = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
        for row in assessment_rows
    )
    locus_cards = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
        for row in locus_rows
    )
    fig_links = "".join(
        f"<li><a href='figures/{html.escape(path.name)}'>{html.escape(path.name)}</a></li>"
        for path in [manhattan] + ([qq] if qq else []) + regional_paths + ld_paths + [evidence_pdf] + list(domain_paths.values())
    )
    html_path.write_text(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Candidate-gene report: {html.escape(trait)}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th, td {{ border: 1px solid #ddd; padding: 6px; vertical-align: top; }}
th {{ background: #f4f6f8; }}
.code {{ font-family: Consolas, monospace; }}
</style>
</head>
<body>
<h1>Candidate-gene report: {html.escape(trait)}</h1>
<p>V3 reports candidate <strong>Priority A/B/C/D</strong> and <strong>Evidence L0-L4</strong> only; it does not claim causal-gene confirmation.</p>
<p>Lead SNPs are <strong>nominal leads</strong> (minimum P within each locus); conditional analysis for multiple independent signals within one locus is not implemented.</p>
<ul>
<li>Data mode: <span class="code">{data_mode}</span></li>
<li>Reference genome: <span class="code">{html.escape(reference_build)}</span></li>
<li>Threshold: <span class="code">{threshold:g}</span></li>
<li>Loci: {len(loci)}; candidate evidence rows: {len(candidates)}</li>
</ul>
<h2>Loci (ranked by nominal lead P)</h2>
<table>
<thead><tr><th>Rank</th><th>Locus</th><th>Lead SNP (pos)</th><th>Lead P</th><th>Cluster type</th><th>LD locus status</th><th>Merged by LD</th></tr></thead>
<tbody>{locus_cards}</tbody>
</table>
<h2>Figures</h2>
<ul>{fig_links}</ul>
<h2>Candidate assessments</h2>
<table>
<thead><tr><th>Gene</th><th>Priority</th><th>Evidence</th><th>TE</th><th>Cloning</th><th>Hormone rescue</th><th>Recommendation</th><th>Downgrade reasons</th></tr></thead>
<tbody>{cards}</tbody>
</table>
</body>
</html>
""", encoding="utf-8")
    return {
        "html": html_path,
        "markdown": md,
        "tables": tables,
        "figures": figures,
        "ld_tables": ld_tables,
    }
