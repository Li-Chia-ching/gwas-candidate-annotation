import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.gates import assess_candidate, detect_te_risk
from gwas_candidate.schemas import HomologyHit


class ExperimentalGateMigrationTest(unittest.TestCase):
    def hit(self):
        return HomologyHit(
            query_gene="GeneA",
            database="ModelPlant",
            hit_id="AT1G00001",
            identity_pct=90.0,
            alignment_length=100,
            query_cover_pct=95.0,
            subject_cover_pct=95.0,
            evalue=1e-20,
            bitscore=200.0,
        )

    def test_te_downgrades_cloning_to_kasp(self):
        self.assertTrue(detect_te_risk(["RVT_2"]))
        assessment = assess_candidate(
            "GeneTE",
            protein_qc_ok=True,
            domains=("RVT_2",),
            homology_hits=(self.hit(),),
        )
        self.assertEqual(assessment.te_status, "detected")
        self.assertEqual(assessment.priority, "C")
        self.assertEqual(assessment.cloning_gate, "forbidden")
        self.assertIn("KASP", assessment.recommendation)
        self.assertIn("te_or_repeat_associated", assessment.downgrade_reasons)

    def test_hormone_rescue_requires_matching_pathway(self):
        mismatch = assess_candidate(
            "GeneA", True, domains=("F-box",), pathway="JA", requested_hormone="GA"
        )
        self.assertEqual(mismatch.hormone_rescue_gate, "forbidden")
        match = assess_candidate(
            "GeneA", True, domains=("GRAS",), pathway="GA", requested_hormone="GA"
        )
        self.assertEqual(match.hormone_rescue_gate, "allowed")

    def test_evidence_and_priority_do_not_overclaim(self):
        weak = assess_candidate("GeneWeak", True)
        self.assertEqual((weak.priority, weak.evidence_level), ("C", "L0"))
        annotated = assess_candidate("GeneA", True, homology_hits=(self.hit(),))
        self.assertEqual((annotated.priority, annotated.evidence_level), ("B", "L2"))
        strong = assess_candidate(
            "GeneA", True, homology_hits=(self.hit(),),
            independent_replication=True, functional_validation=True,
        )
        self.assertEqual((strong.priority, strong.evidence_level), ("A", "L4"))
        failed_qc = assess_candidate("GeneBad", False)
        self.assertEqual(failed_qc.priority, "D")
        self.assertEqual(failed_qc.cloning_gate, "forbidden")


if __name__ == "__main__":
    unittest.main()
