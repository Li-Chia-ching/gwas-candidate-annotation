"""Multi-database BLASTp and homology-candidate adapters.

BLAST top hits are homology candidates only. They are not ortholog calls and are not
true homoeolog copy counts.
"""

import csv
import subprocess
from dataclasses import replace
from pathlib import Path

from .schemas import HomologyHit


OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qcovs qlen slen"
)


def read_db_config(path):
    path = Path(path)
    databases = []
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected name<TAB>protein_fasta")
            fasta = Path(fields[1])
            if not fasta.is_absolute():
                fasta = path.parent / fasta
            databases.append((fields[0], fasta))
    return databases


def run_multi_db_blastp(
    query_fasta,
    db_config,
    outdir,
    blastp="blastp",
    makeblastdb="makeblastdb",
    evalue="1e-5",
    max_target=20,
    threads=1,
):
    outdir = Path(outdir)
    blast_out = outdir / "blast_out"
    db_root = outdir / "db"
    blast_out.mkdir(parents=True, exist_ok=True)
    db_root.mkdir(parents=True, exist_ok=True)
    outputs = []
    for name, fasta in read_db_config(db_config):
        if not fasta.exists():
            raise FileNotFoundError(f"database FASTA missing for {name}: {fasta}")
        db_dir = db_root / name
        db_dir.mkdir(parents=True, exist_ok=True)
        prefix = db_dir / name
        psq = db_dir / f"{name}.psq"
        if not psq.exists() or psq.stat().st_mtime < fasta.stat().st_mtime:
            subprocess.run(
                [makeblastdb, "-in", str(fasta), "-dbtype", "prot", "-out", str(prefix), "-parse_seqids"],
                check=True,
            )
        output = blast_out / f"{name}.blast.tsv"
        subprocess.run(
            [
                blastp, "-query", str(query_fasta), "-db", str(prefix), "-out", str(output),
                "-outfmt", OUTFMT, "-evalue", str(evalue),
                "-max_target_seqs", str(max_target), "-num_threads", str(threads),
            ],
            check=True,
        )
        outputs.append(output)
    return outputs


def parse_blast_hits(path, database):
    hits = []
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 15:
                raise ValueError(f"{path}:{line_no}: expected 15 BLAST columns")
            query_gene = fields[0].split("|")[0]
            subject_length = float(fields[14])
            subject_cover = (abs(int(fields[9]) - int(fields[8])) + 1) / subject_length * 100
            hits.append(HomologyHit(
                query_gene=query_gene,
                database=database,
                hit_id=fields[1],
                identity_pct=float(fields[2]),
                alignment_length=int(fields[3]),
                query_cover_pct=float(fields[12]),
                subject_cover_pct=subject_cover,
                evalue=float(fields[10]),
                bitscore=float(fields[11]),
            ))
    return hits


def top_homology_candidates(hits, top=5):
    grouped = {}
    for hit in sorted(hits, key=lambda h: (h.query_gene, h.database, h.evalue, -h.bitscore)):
        key = (hit.query_gene, hit.database)
        seen = {h.hit_id for h in grouped.get(key, [])}
        if hit.hit_id not in seen:
            grouped.setdefault(key, []).append(hit)
    return [hit for key in sorted(grouped) for hit in grouped[key][:top]]


def self_blast_candidates(hits, self_database):
    return [
        replace(hit, evidence_label="homoeolog_candidate")
        for hit in hits if hit.database == self_database
    ]


def write_homology_tables(hits, outdir, self_database=None, top=5):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    top_hits = top_homology_candidates(hits, top=top)
    summary = outdir / "homology_candidates.tsv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_gene", "database", "hit_id", "identity_pct", "alignment_length",
            "query_cover_pct", "subject_cover_pct", "evalue", "bitscore", "evidence_label"
        ])
        for hit in top_hits:
            writer.writerow([
                hit.query_gene, hit.database, hit.hit_id, hit.identity_pct,
                hit.alignment_length, hit.query_cover_pct,
                "" if hit.subject_cover_pct is None else hit.subject_cover_pct,
                hit.evalue, hit.bitscore, hit.evidence_label,
            ])
    self_path = None
    if self_database:
        self_hits = self_blast_candidates(top_hits, self_database)
        self_path = outdir / "homoeolog_candidates.tsv"
        with self_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "query_gene", "candidate_gene", "identity_pct", "alignment_length",
                "query_cover_pct", "subject_cover_pct", "evalue", "evidence_label"
            ])
            for hit in self_hits:
                writer.writerow([
                    hit.query_gene, hit.hit_id, hit.identity_pct, hit.alignment_length,
                    hit.query_cover_pct, hit.subject_cover_pct, hit.evalue, hit.evidence_label,
                ])
    return summary, self_path
