"""gwas-candidate-annotation V3 adapter package."""

from .association import build_dual_model_clusters, read_model_signals
from .candidates import candidate_genes_for_loci, candidate_genes_for_locus, load_gff_genes
from .gates import assess_candidate, detect_te_risk
from .homology import parse_blast_hits, self_blast_candidates, top_homology_candidates
from .ld import ld_matrix_for_region
from .pipeline import run_candidate_pipeline
from .reporting import generate_report
from .sequence import ProteinQCError, extract_cds_protein, extract_proteome

__all__ = [
    "ProteinQCError",
    "assess_candidate",
    "build_dual_model_clusters",
    "candidate_genes_for_loci",
    "candidate_genes_for_locus",
    "detect_te_risk",
    "extract_cds_protein",
    "extract_proteome",
    "generate_report",
    "ld_matrix_for_region",
    "load_gff_genes",
    "parse_blast_hits",
    "read_model_signals",
    "run_candidate_pipeline",
    "self_blast_candidates",
    "top_homology_candidates",
]
