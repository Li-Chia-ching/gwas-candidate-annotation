"""End-to-end candidate annotation pipeline over the V3 adapter interfaces."""

import csv
from dataclasses import replace
from pathlib import Path

from .association import build_dual_model_clusters, read_model_signals
from .candidates import candidate_genes_for_loci, load_gff_genes
from .gates import assess_candidate
from .homology import parse_blast_hits, run_multi_db_blastp, top_homology_candidates
from .provenance import write_provenance
from .reporting import generate_report
from .sequence import extract_cds_protein, write_sequence_outputs


def _resolve(base, value):
    path = Path(value)
    return path if path.is_absolute() else Path(base) / path


def read_domains(path):
    domains = {}
    if not path:
        return domains
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "name", "start", "end"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"domain table requires {sorted(required)}")
        for row in reader:
            domains.setdefault(row["gene_id"], []).append({
                "name": row["name"], "start": int(row["start"]), "end": int(row["end"])
            })
    return domains


def run_candidate_pipeline(config, workdir):
    """Run association adapters through candidate evidence reports."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    base_value = config.get("base_dir")
    base = _resolve(workdir, base_value) if base_value else workdir
    trait = config["trait"]
    data_mode = config["data_mode"]
    reference_build = config["reference_build"]
    commands = []

    model_a = config["model_a"]
    model_b = config["model_b"]
    signals_a_path = _resolve(base, model_a["path"])
    signals_b_path = _resolve(base, model_b["path"])
    signals_a = read_model_signals(
        signals_a_path, model_a["name"], trait, float(model_a["p_threshold"]), data_mode
    )
    signals_b = read_model_signals(
        signals_b_path, model_b["name"], trait, float(model_b["p_threshold"]), data_mode
    )
    commands.append(
        f"read_model_signals({signals_a_path.name}, {signals_b_path.name}, "
        f"thresholds={model_a['p_threshold']},{model_b['p_threshold']})"
    )
    loci = build_dual_model_clusters(
        signals_a, signals_b, model_a["name"], model_b["name"], trait,
        int(config.get("gap", 100000)), reference_build, data_mode,
    )

    gff = _resolve(base, config["gff"])
    genome = _resolve(base, config["genome"])
    genes = load_gff_genes(gff)
    candidates = candidate_genes_for_loci(loci, genes, int(config.get("flank", 100000)))
    gene_ids = sorted({c.gene_id for c in candidates})
    proteins = extract_cds_protein(gff, genome, gene_ids, strict=True)
    sequence_dir = workdir / "sequences"
    write_sequence_outputs(proteins, sequence_dir)
    commands.append("extract_cds_protein(strict=True)")

    hits = []
    db_config = config.get("db_config")
    if db_config:
        query = sequence_dir / "candidate_proteins.fasta"
        blast_paths = run_multi_db_blastp(
            query,
            _resolve(base, db_config),
            workdir / "homology",
            evalue=str(config.get("evalue", "1e-5")),
            max_target=int(config.get("max_target", 20)),
            threads=1,
        )
        commands.append("run_multi_db_blastp(local BLAST+)")
        for path in blast_paths:
            hits.extend(parse_blast_hits(path, path.name.replace(".blast.tsv", "")))
    hits = top_homology_candidates(hits, top=int(config.get("top_hits", 5)))
    self_db = config.get("self_db")
    if self_db:
        hits = [
            replace(hit, evidence_label="homoeolog_candidate") if hit.database == self_db else hit
            for hit in hits
        ]

    domains = read_domains(_resolve(base, config["domains"]) if config.get("domains") else None)
    pathways = config.get("pathways", {})
    requested_hormone = config.get("requested_hormone")
    assessments = {}
    for gene_id in gene_ids:
        gene_hits = tuple(h for h in hits if h.query_gene == gene_id)
        assessments[gene_id] = assess_candidate(
            gene_id,
            protein_qc_ok=True,
            domains=domains.get(gene_id, ()),
            homology_hits=gene_hits,
            pathway=pathways.get(gene_id),
            requested_hormone=requested_hormone,
            ld_support=data_mode == "full" and bool(config.get("vcf")),
        )

    report = generate_report(
        trait=trait,
        association_path=signals_a_path,
        loci=loci,
        candidates=candidates,
        proteins=proteins,
        homology_hits=hits,
        assessments=assessments,
        domains=domains,
        outdir=workdir / "report",
        threshold=float(model_a["p_threshold"]),
        data_mode=data_mode,
        reference_build=reference_build,
        vcf_path=_resolve(base, config["vcf"]) if config.get("vcf") else None,
        region_flank=int(config.get("region_flank", config.get("flank", 100000))),
    )
    inputs = [signals_a_path, signals_b_path, gff, genome]
    for key in ("vcf", "db_config", "domains"):
        if config.get(key):
            inputs.append(_resolve(base, config[key]))
    manifest = write_provenance(
        workdir, config.get("scenario", trait), config, inputs, commands,
        int(config.get("random_seed", 0)), reference_build,
    )
    return {
        "loci": loci,
        "candidates": candidates,
        "proteins": proteins,
        "homology_hits": hits,
        "assessments": assessments,
        "report": report,
        "manifest": manifest,
    }
