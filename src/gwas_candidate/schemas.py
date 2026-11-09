"""Shared schemas for gwas-candidate-annotation V3."""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class VariantRecord:
    trait: str
    model: str
    chrom: str
    pos: int
    snp: str
    ref: str
    alt: str
    p: float
    effect: float
    data_mode: str
    annotations: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class LocusRecord:
    trait: str
    chrom: str
    start: int
    end: int
    lead: VariantRecord
    lead_model: str
    support_models: tuple[str, ...]
    variants: tuple[VariantRecord, ...]
    reference_build: str
    data_mode: str
    cluster_type: str = "distance_cluster"
    ld_locus_status: str = "not_evaluated"


@dataclass(frozen=True)
class GeneRecord:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class CandidateGeneRecord:
    trait: str
    locus_id: str
    chrom: str
    lead_pos: int
    gene_id: str
    gene_start: int
    gene_end: int
    strand: str
    distance_to_lead_bp: int
    overlaps_lead: bool
    evidence_source: str


@dataclass(frozen=True)
class ProteinRecord:
    gene_id: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    cds_sequence: str
    protein_sequence: str
    qc_status: str
    qc_issues: tuple[str, ...] = ()


@dataclass(frozen=True)
class HomologyHit:
    query_gene: str
    database: str
    hit_id: str
    identity_pct: float
    alignment_length: int
    query_cover_pct: float
    subject_cover_pct: float | None
    evalue: float
    bitscore: float
    evidence_label: str = "homology_candidate"


@dataclass(frozen=True)
class CandidateAssessment:
    gene_id: str
    priority: str
    evidence_level: str
    te_status: str
    cloning_gate: str
    hormone_rescue_gate: str
    recommendation: str
    downgrade_reasons: tuple[str, ...] = ()
