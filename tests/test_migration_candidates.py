import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.association import build_dual_model_clusters, read_model_signals
from gwas_candidate.candidates import candidate_genes_for_locus, load_gff_genes


class CandidateAdapterMigrationTest(unittest.TestCase):
    def test_dual_model_adapter_preserves_primary_lead_and_cluster_contract(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            files = []
            for name, p1, p2 in (("a", "1e-6", "1e-9"), ("b", "1e-8", "1e-7")):
                path = tmp / f"{name}.csv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect", "snpEff_gene"])
                    writer.writerow(["s1", "chr1", "100", p1, "A", "G", "0.2", "GeneA"])
                    writer.writerow(["s2", "chr1", "105", p2, "C", "T", "0.4", ""])
                files.append(path)
            sig_a = read_model_signals(files[0], "MLM", "T", 1e-5, "full")
            sig_b = read_model_signals(files[1], "FarmCPU", "T", 1e-5, "full")
            loci = build_dual_model_clusters(
                sig_a, sig_b, "MLM", "FarmCPU", "T", 10, "buildA", "full"
            )
            self.assertEqual(len(loci), 1)
            self.assertEqual(loci[0].lead.pos, 105)
            self.assertEqual(loci[0].cluster_type, "distance_cluster")
            self.assertEqual(loci[0].ld_locus_status, "not_evaluated")
            self.assertEqual(loci[0].support_models, ("MLM", "FarmCPU"))

    def test_gff_adapter_preserves_window_direct_and_gene_body_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff = tmp / "ref.gff"
            gff.write_text(
                "chr1\ttest\tgene\t80\t120\t.\t+\t.\tID=GeneA\n"
                "chr1\ttest\tgene\t150\t180\t.\t-\t.\tID=GeneB\n"
                "chr1\ttest\tgene\t500\t600\t.\t+\t.\tID=GeneC\n",
                encoding="utf-8",
            )
            genes = load_gff_genes(gff)
            self.assertEqual(genes["GeneB"].strand, "-")
            path = tmp / "a.csv"
            path.write_text(
                "SNP,CHROM,POS,P.value,REF,ALT,Effect,snpEff_gene\n"
                "s1,chr1,100,1e-8,A,G,0.2,GeneA\n",
                encoding="utf-8",
            )
            signals = read_model_signals(path, "MLM", "T", 1e-5, "full")
            locus = build_dual_model_clusters(
                signals, signals, "MLM", "MLM", "T", 10, "buildA", "full"
            )[0]
            records = candidate_genes_for_locus(locus, genes, flank=100)
            by_source = {(r.gene_id, r.evidence_source): r for r in records}
            self.assertIn(("GeneA", "direct_annotation"), by_source)
            self.assertIn(("GeneA", "gene_body_overlap"), by_source)
            self.assertIn(("GeneB", "fixed_window"), by_source)
            self.assertTrue(by_source[("GeneA", "gene_body_overlap")].overlaps_lead)
            self.assertEqual(by_source[("GeneB", "fixed_window")].distance_to_lead_bp, 50)
            self.assertNotIn("GeneC", {r.gene_id for r in records})


if __name__ == "__main__":
    unittest.main()
