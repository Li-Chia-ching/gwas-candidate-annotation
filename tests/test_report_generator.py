import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.gates import assess_candidate
from gwas_candidate.reporting import generate_report
from gwas_candidate.schemas import (
    CandidateGeneRecord,
    HomologyHit,
    LocusRecord,
    ProteinRecord,
    VariantRecord,
)


class CandidateReportTest(unittest.TestCase):
    def make_bundle(self, tmp):
        assoc = tmp / "assoc.tsv"
        assoc.write_text(
            "CHR\tBP\tP\tBETA\tN\n"
            "1\t100\t1e-8\t0.2\t4\n"
            "1\t110\t2e-8\t0.3\t4\n"
            "1\t500\t0.2\t0.1\t4\n",
            encoding="utf-8",
        )
        vcf = tmp / "ld.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
            "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t1/1\n"
            "chr1\t110\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t1/1\n",
            encoding="utf-8",
        )
        variant = VariantRecord(
            trait="height", model="MLM", chrom="chr1", pos=100, snp="s1",
            ref="A", alt="G", p=1e-8, effect=0.2, data_mode="full",
            annotations={"hgvs_p": "p.Ala2Val", "gene_ann": "GeneA"},
        )
        locus = LocusRecord(
            trait="height", chrom="chr1", start=100, end=110, lead=variant,
            lead_model="MLM", support_models=("MLM", "FarmCPU"), variants=(variant,),
            reference_build="buildA", data_mode="full",
        )
        candidates = [
            CandidateGeneRecord(
                trait="height", locus_id="height:chr1:100-110", chrom="chr1",
                lead_pos=100, gene_id="GeneA", gene_start=80, gene_end=120,
                strand="+", distance_to_lead_bp=0, overlaps_lead=True,
                evidence_source="direct_annotation",
            ),
            CandidateGeneRecord(
                trait="height", locus_id="height:chr1:100-110", chrom="chr1",
                lead_pos=100, gene_id="GeneA", gene_start=80, gene_end=120,
                strand="+", distance_to_lead_bp=0, overlaps_lead=True,
                evidence_source="gene_body_overlap",
            ),
        ]
        proteins = [ProteinRecord(
            gene_id="GeneA", transcript_id="TxA", chrom="chr1", start=80,
            end=120, strand="+", cds_sequence="ATGGCTTAA",
            protein_sequence="MA*", qc_status="OK",
        )]
        hits = [HomologyHit(
            query_gene="GeneA", database="ModelPlant", hit_id="AT1G00001",
            identity_pct=95.0, alignment_length=2, query_cover_pct=100.0,
            subject_cover_pct=100.0, evalue=1e-10, bitscore=50.0,
        )]
        assessments = {"GeneA": assess_candidate("GeneA", True, domains=("DomainA",), homology_hits=hits)}
        domains = {"GeneA": [{"name": "DomainA", "start": 1, "end": 2}]}
        return assoc, vcf, [locus], candidates, proteins, hits, assessments, domains

    def test_full_report_outputs_html_markdown_tables_and_pdf_figures(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            bundle = self.make_bundle(tmp)
            result = generate_report(
                "height", bundle[0], bundle[2], bundle[3], bundle[4], bundle[5],
                bundle[6], bundle[7], tmp / "report", threshold=1e-5,
                data_mode="full", reference_build="buildA", vcf_path=bundle[1],
                region_flank=100,
            )
            expected = [
                result["html"], result["markdown"],
                result["tables"] / "loci.tsv",
                result["tables"] / "candidate_genes.tsv",
                result["tables"] / "homology_candidates.tsv",
                result["tables"] / "candidate_assessments.tsv",
                result["figures"] / "height_manhattan.pdf",
                result["figures"] / "height_qq.pdf",
                result["figures"] / "locus_01_regional_association.pdf",
                result["figures"] / "locus_01_ld_heatmap.pdf",
                result["figures"] / "candidate_evidence_matrix.pdf",
                result["figures"] / "domain_map_GeneA.pdf",
            ]
            for path in expected:
                self.assertTrue(path.exists(), path)
            for pdf in [p for p in expected if p.suffix == ".pdf"]:
                self.assertEqual(pdf.read_bytes()[:4], b"%PDF")
            html = result["html"].read_text(encoding="utf-8")
            self.assertIn("Priority", html)
            self.assertIn("Evidence", html)
            self.assertIn("homology_candidate", (result["tables"] / "homology_candidates.tsv").read_text(encoding="utf-8"))

    def test_truncated_report_omits_qq(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            bundle = self.make_bundle(tmp)
            result = generate_report(
                "height", bundle[0], bundle[2], bundle[3], bundle[4], bundle[5],
                bundle[6], bundle[7], tmp / "report", threshold=1e-5,
                data_mode="truncated", reference_build="buildA", vcf_path=bundle[1],
                region_flank=100,
            )
            self.assertFalse((result["figures"] / "height_qq.pdf").exists())
            self.assertIn("not evaluated for truncated input", result["markdown"].read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
