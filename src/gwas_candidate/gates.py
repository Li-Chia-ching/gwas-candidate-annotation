"""Experimental gating for candidate genes."""

import re

from .schemas import CandidateAssessment, HomologyHit


TE_PATTERN = re.compile(r"(RVT_|rve|Retrotrans_gag|reverse transcriptase|transpos)", re.IGNORECASE)


def detect_te_risk(domains=(), homology_descriptions=()):
    text = "\n".join([str(x) for x in domains] + [str(x) for x in homology_descriptions])
    return bool(TE_PATTERN.search(text))


def assess_candidate(
    gene_id,
    protein_qc_ok,
    domains=(),
    homology_hits=(),
    pathway=None,
    requested_hormone=None,
    ld_support=False,
    independent_replication=False,
    functional_validation=False,
):
    """Assign Priority A-D and Evidence L0-L4 without claiming causal confirmation."""
    domains = tuple(domains)
    hits = tuple(homology_hits)
    te = detect_te_risk(domains)
    reasons = []
    if not protein_qc_ok:
        reasons.append("protein_qc_failed")
    if te:
        reasons.append("te_or_repeat_associated")

    if functional_validation:
        evidence_level = "L4"
    elif independent_replication:
        evidence_level = "L3"
    elif domains or hits:
        evidence_level = "L2"
    elif ld_support:
        evidence_level = "L1"
    else:
        evidence_level = "L0"

    if not protein_qc_ok:
        priority = "D"
        cloning_gate = "forbidden"
    elif te:
        priority = "C"
        cloning_gate = "forbidden"
    elif functional_validation and independent_replication:
        priority = "A"
        cloning_gate = "allowed"
    elif independent_replication or domains or hits:
        priority = "B"
        cloning_gate = "allowed"
    else:
        priority = "C"
        cloning_gate = "needs_evidence"

    if requested_hormone is None:
        hormone_gate = "not_requested"
    elif pathway and pathway.lower() == requested_hormone.lower():
        hormone_gate = "allowed"
    else:
        hormone_gate = "forbidden"
        reasons.append("hormone_pathway_not_established")

    if not protein_qc_ok:
        recommendation = "repair gene model and repeat CDS/protein QC before experimental design"
    elif te:
        recommendation = "KASP marker validation only; cloning is downgraded"
    elif hormone_gate == "allowed":
        recommendation = "clone with copy-aware controls; pathway-matched hormone rescue is allowed"
    elif cloning_gate == "allowed":
        recommendation = "clone after marker validation; add independent replication evidence"
    else:
        recommendation = "collect functional evidence before cloning"

    return CandidateAssessment(
        gene_id=gene_id,
        priority=priority,
        evidence_level=evidence_level,
        te_status="detected" if te else "not_detected",
        cloning_gate=cloning_gate,
        hormone_rescue_gate=hormone_gate,
        recommendation=recommendation,
        downgrade_reasons=tuple(reasons),
    )
