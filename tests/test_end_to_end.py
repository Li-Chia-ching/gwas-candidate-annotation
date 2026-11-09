import json
import shutil
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.association import build_dual_model_clusters, read_model_signals
from gwas_candidate.pipeline import run_candidate_pipeline
from gwas_candidate.provenance import write_provenance
from gwas_candidate.sequence import revcomp


E2E_ROOT = ROOT / "artifacts" / "e2e"
CODONS = [
    "ATG", "GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAG",
    "TTG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT",
]
CODING = "".join(CODONS) + "TAA"
PROTEIN = "MACDEFGHIKLNPQRSTVWY*"


class EndToEndScenariosTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if E2E_ROOT.exists():
            shutil.rmtree(E2E_ROOT)
        E2E_ROOT.mkdir(parents=True)
        cls.summary = []

    @classmethod
    def tearDownClass(cls):
        lines = [
            "# V3 constructed-data E2E summary",
            "",
            "All scenarios use constructed data, random_seed=20261109, reference_genome_version=construct-v1.",
            "",
            "| Scenario | Status | Evidence |",
            "| --- | --- | --- |",
        ]
        lines += cls.summary
        lines += [
            "",
            "BLAST labels are limited to homology_candidate and homoeolog_candidate; no ortholog or true copy-count claim is emitted.",
        ]
        (E2E_ROOT / "E2E_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def make_scenario(self, name, negative=False, zero_blast=False, polyploid=False, truncated=False):
        scenario = E2E_ROOT / name
        inputs = scenario / "inputs"
        inputs.mkdir(parents=True)
        (inputs / "model_a.csv").write_text(
            "SNP,CHROM,POS,P.value,REF,ALT,Effect,snpEff_gene,hgvs_p\n"
            "s1,chr1,10,1e-8,A,G,0.2,GeneA,p.Ala2Val\n"
            "s2,chr1,20,2e-8,C,T,0.3,GeneA,\n"
            "s3,chr1,500,1e-3,G,A,0.1,,\n",
            encoding="utf-8",
        )
        (inputs / "model_b.csv").write_text(
            "SNP,CHROM,POS,P.value,REF,ALT,Effect\n"
            "s1,chr1,10,2e-8,A,G,0.4\n"
            "s2,chr1,20,3e-8,C,T,0.5\n"
            "s3,chr1,500,4e-8,G,A,0.2\n",
            encoding="utf-8",
        )
        end = len(CODING)
        if polyploid:
            coding_b = CODING
            gene_b_start = end + 17
            gene_b_end = gene_b_start + len(coding_b) - 1
            genome_seq = CODING + "N" * 16 + coding_b + "N" * 20
            gff_text = (
                f"chr1\tconstruct\tgene\t1\t{end}\t.\t+\t.\tID=GeneA\n"
                f"chr1\tconstruct\tmRNA\t1\t{end}\t.\t+\t.\tID=TxA;Parent=GeneA\n"
                f"chr1\tconstruct\tCDS\t1\t{end}\t.\t+\t0\tParent=TxA\n"
                f"chr1\tconstruct\tgene\t{gene_b_start}\t{gene_b_end}\t.\t+\t.\tID=GeneB\n"
                f"chr1\tconstruct\tmRNA\t{gene_b_start}\t{gene_b_end}\t.\t+\t.\tID=TxB;Parent=GeneB\n"
                f"chr1\tconstruct\tCDS\t{gene_b_start}\t{gene_b_end}\t.\t+\t0\tParent=TxB\n"
            )
            proteome = ">GeneA\n" + PROTEIN + "\n>GeneB\n" + PROTEIN + "\n"
            (inputs / "proteome.fa").write_text(proteome, encoding="utf-8")
            db_text = "Self\tproteome.fa\n"
        else:
            genome_seq = (revcomp(CODING) if negative else CODING) + "N" * 37
            strand = "-" if negative else "+"
            gff_text = (
                f"chr1\tconstruct\tgene\t1\t{end}\t.\t{strand}\t.\tID=GeneA\n"
                f"chr1\tconstruct\tmRNA\t1\t{end}\t.\t{strand}\t.\tID=TxA;Parent=GeneA\n"
                f"chr1\tconstruct\tCDS\t1\t{end}\t.\t{strand}\t0\tParent=TxA\n"
            )
            db_protein = "K" * 20 if zero_blast else PROTEIN
            (inputs / "model_plant.fa").write_text(">ModelProtein\n" + db_protein + "\n", encoding="utf-8")
            db_text = "ModelPlant\tmodel_plant.fa\n"
        (inputs / "genome.fa").write_text(">chr1\n" + genome_seq + "\n", encoding="utf-8")
        (inputs / "ref.gff").write_text(gff_text, encoding="utf-8")
        (inputs / "dbs.tsv").write_text(db_text, encoding="utf-8")
        gt_values = (
            ["0/0/0/0", "0/0/0/1", "1/1/1/1", "1/1/1/1"]
            if polyploid else ["0/0", "0/1", "1/1", "1/1"]
        )
        vcf_lines = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
            "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t" + "\t".join(gt_values),
            "chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t" + "\t".join(gt_values),
        ]
        (inputs / "variants.vcf").write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")
        (inputs / "domains.tsv").write_text(
            "" if zero_blast else "gene_id\tname\tstart\tend\nGeneA\tAla_rich\t2\t8\n",
            encoding="utf-8",
        )
        config = {
            "scenario": name,
            "base_dir": "inputs",
            "trait": "height",
            "data_mode": "truncated" if truncated else "full",
            "reference_build": "construct-v1",
            "random_seed": 20261109,
            "model_a": {"name": "MLM", "path": "model_a.csv", "p_threshold": 1e-5},
            "model_b": {"name": "FarmCPU", "path": "model_b.csv", "p_threshold": 1e-5},
            "gff": "ref.gff",
            "genome": "genome.fa",
            "db_config": "dbs.tsv",
            "vcf": "variants.vcf",
            "gap": 20,
            "flank": 50,
            "region_flank": 30,
            "evalue": "1e-50" if zero_blast else "1e-5",
        }
        if not zero_blast:
            config["domains"] = "domains.tsv"
        if polyploid:
            config["self_db"] = "Self"
            config["flank"] = 100
        return scenario, config

    def record(self, scenario, status, evidence):
        self.summary.append(f"| {scenario} | {status} | {evidence} |")

    def test_01_full_data_end_to_end(self):
        scenario, config = self.make_scenario("full")
        result = run_candidate_pipeline(config, scenario)
        self.assertTrue((result["report"]["figures"] / "height_qq.pdf").exists())
        self.assertTrue((result["report"]["figures"] / "locus_01_ld_heatmap.pdf").exists())
        self.assertTrue(result["manifest"].exists())
        manifest_text = result["manifest"].read_text(encoding="utf-8")
        report_text = result["report"]["markdown"].read_text(encoding="utf-8")
        self.assertNotIn(str(scenario), manifest_text)
        self.assertNotIn(str(scenario), report_text)
        self.assertNotRegex(manifest_text, r"C:[/\\]+Users")
        self.assertNotRegex(report_text, r"C:[/\\]+Users")
        self.assertEqual(result["assessments"]["GeneA"].evidence_level, "L2")
        self.assertTrue(result["homology_hits"])
        (scenario / "scenario.log").write_text("PASS: full mode generated QQ, LD heatmap, BLAST hit, and report.\n", encoding="utf-8")
        self.record("full", "PASS", "QQ + LD heatmap + BLAST homology candidate + report")

    def test_02_truncated_data_end_to_end(self):
        scenario, config = self.make_scenario("truncated", truncated=True)
        result = run_candidate_pipeline(config, scenario)
        self.assertFalse((result["report"]["figures"] / "height_qq.pdf").exists())
        self.assertTrue((result["report"]["figures"] / "height_manhattan.pdf").exists())
        self.assertIn("not evaluated for truncated", (result["report"]["markdown"]).read_text(encoding="utf-8"))
        (scenario / "scenario.log").write_text("PASS: truncated mode generated Manhattan/report and omitted QQ/lambda.\n", encoding="utf-8")
        self.record("truncated", "PASS", "QQ/lambda disabled; report generated")

    def test_03_multiallelic_data_blocks(self):
        scenario = E2E_ROOT / "multiallelic"
        inputs = scenario / "inputs"
        inputs.mkdir(parents=True)
        (inputs / "a.csv").write_text(
            "SNP,CHROM,POS,P.value,REF,ALT,Effect\ns1,chr1,10,1e-8,A,G,0.2\n",
            encoding="utf-8",
        )
        (inputs / "b.csv").write_text(
            "SNP,CHROM,POS,P.value,REF,ALT,Effect\ns1,chr1,10,1e-8,A,T,0.2\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "allele mismatch"):
            sig_a = read_model_signals(inputs / "a.csv", "MLM", "height", 1e-5, "full")
            sig_b = read_model_signals(inputs / "b.csv", "FarmCPU", "height", 1e-5, "full")
            build_dual_model_clusters(sig_a, sig_b, "MLM", "FarmCPU", "height", 20, "construct-v1", "full")
        config = {"scenario": "multiallelic", "data_mode": "full"}
        write_provenance(
            scenario, "multiallelic", config, [inputs / "a.csv", inputs / "b.csv"],
            ["read_model_signals + build_dual_model_clusters"], 20261109, "construct-v1",
        )
        (scenario / "scenario.log").write_text("PASS: multiallelic allele conflict blocked.\n", encoding="utf-8")
        self.record("multiallelic", "PASS", "expected hard failure: allele mismatch")

    def test_04_negative_strand_gene_end_to_end(self):
        scenario, config = self.make_scenario("negative_strand", negative=True)
        result = run_candidate_pipeline(config, scenario)
        self.assertEqual(result["proteins"][0].protein_sequence, PROTEIN)
        self.assertEqual(result["proteins"][0].qc_status, "OK")
        self.assertTrue(result["report"]["html"].exists())
        (scenario / "scenario.log").write_text("PASS: negative-strand CDS translated correctly and reported.\n", encoding="utf-8")
        self.record("negative_strand", "PASS", "reverse-complement CDS + report")

    def test_05_zero_blast_hit_end_to_end(self):
        scenario, config = self.make_scenario("zero_blast", zero_blast=True)
        result = run_candidate_pipeline(config, scenario)
        self.assertEqual(result["homology_hits"], [])
        self.assertEqual(result["assessments"]["GeneA"].evidence_level, "L1")
        text = (result["report"]["tables"] / "homology_candidates.tsv").read_text(encoding="utf-8")
        self.assertNotIn("ortholog", text.lower())
        self.assertTrue(result["manifest"].exists())
        (scenario / "scenario.log").write_text("PASS: zero BLAST hit retained with explicit L1 LD-only evidence.\n", encoding="utf-8")
        self.record("zero_blast", "PASS", "zero-hit table + L1/C assessment")

    def test_06_polyploid_self_blast_end_to_end(self):
        scenario, config = self.make_scenario("polyploid", polyploid=True)
        result = run_candidate_pipeline(config, scenario)
        self.assertEqual(len(result["proteins"]), 2)
        labels = {hit.evidence_label for hit in result["homology_hits"] if hit.database == "Self"}
        self.assertEqual(labels, {"homoeolog_candidate"})
        self.assertTrue((result["report"]["figures"] / "locus_01_ld_heatmap.pdf").exists())
        text = (result["report"]["tables"] / "homology_candidates.tsv").read_text(encoding="utf-8")
        self.assertNotIn("copy_number", text.lower())
        (scenario / "scenario.log").write_text("PASS: tetraploid dosage LD and self-BLAST homoeolog candidates generated.\n", encoding="utf-8")
        self.record("polyploid", "PASS", "tetraploid LD + self-BLAST homoeolog candidates")


if __name__ == "__main__":
    unittest.main()
