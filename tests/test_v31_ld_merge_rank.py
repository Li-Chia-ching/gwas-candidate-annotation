import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.association import (
    build_dual_model_clusters,
    merge_clusters_by_ld,
    rank_loci,
    read_model_signals,
)
from gwas_candidate.pipeline import run_candidate_pipeline
from gwas_candidate.schemas import LocusRecord, VariantRecord


CODONS = [
    "ATG", "GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAG",
    "TTG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT",
]
CODING = "".join(CODONS) + "TAA"


def write_signals(path, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect"])
        writer.writerows(rows)


def write_vcf(path, records, chrom="chr1"):
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    ]
    for pos, ref, alt, gts in records:
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_loci(tmp, rows_a, rows_b, gap=20):
    a = tmp / "a.csv"
    b = tmp / "b.csv"
    write_signals(a, rows_a)
    write_signals(b, rows_b)
    signals_a = read_model_signals(a, "MLM", "height", 1e-5, "full")
    signals_b = read_model_signals(b, "FarmCPU", "height", 1e-5, "full")
    return build_dual_model_clusters(
        signals_a, signals_b, "MLM", "FarmCPU", "height", gap, "buildA", "full"
    )


def make_locus(chrom, pos, p):
    variant = VariantRecord(
        trait="height", model="MLM", chrom=chrom, pos=pos, snp=f"s{pos}",
        ref="A", alt="G", p=p, effect=0.2, data_mode="full",
    )
    return LocusRecord(
        trait="height", chrom=chrom, start=pos, end=pos, lead=variant,
        lead_model="MLM", support_models=("MLM", "FarmCPU"), variants=(variant,),
        reference_build="buildA", data_mode="full",
    )


ROWS_A = [
    ["s1", "chr1", "10", "1e-8", "A", "G", "0.2"],
    ["s2", "chr1", "80", "2e-8", "C", "T", "0.3"],
]
ROWS_B = [
    ["s1", "chr1", "10", "3e-8", "A", "G", "0.4"],
    ["s2", "chr1", "80", "4e-8", "C", "T", "0.5"],
]


class LdMergeTest(unittest.TestCase):
    def test_high_ld_merges_adjacent_clusters(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            loci = build_loci(tmp, ROWS_A, ROWS_B)
            self.assertEqual(len(loci), 2)
            vcf = tmp / "ld.vcf"
            write_vcf(vcf, [
                (10, "A", "G", ["0/0", "0/1", "1/1", "1/1"]),
                (80, "C", "T", ["0/0", "0/1", "1/1", "1/1"]),
            ])
            merged = merge_clusters_by_ld(loci, vcf, 0.5)
            self.assertEqual(len(merged), 1)
            locus = merged[0]
            self.assertTrue(locus.merged_by_ld)
            self.assertEqual(locus.merged_cluster_count, 2)
            self.assertEqual(locus.cluster_type, "ld_merged_cluster")
            self.assertEqual(locus.ld_locus_status, "merged_by_ld")
            self.assertEqual(locus.start, 10)
            self.assertEqual(locus.end, 80)
            self.assertEqual(locus.lead.snp, "s1")  # nominal lead = min P
            self.assertEqual(len(locus.variants), 2)
            self.assertEqual(locus.support_models, ("MLM", "FarmCPU"))

    def test_low_ld_keeps_clusters_separate(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            loci = build_loci(tmp, ROWS_A, ROWS_B)
            vcf = tmp / "ld.vcf"
            write_vcf(vcf, [
                (10, "A", "G", ["0/0", "0/0", "1/1", "1/1"]),
                (80, "C", "T", ["1/1", "0/0", "1/1", "0/0"]),
            ])
            result = merge_clusters_by_ld(loci, vcf, 0.5)
            self.assertEqual(len(result), 2)
            for locus in result:
                self.assertFalse(locus.merged_by_ld)
                self.assertEqual(locus.cluster_type, "distance_cluster")
                self.assertEqual(locus.ld_locus_status, "ld_evaluated_not_merged")

    def test_no_vcf_degrades_to_not_evaluated(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            loci = build_loci(tmp, ROWS_A, ROWS_B)
            result = merge_clusters_by_ld(loci, None, 0.5)
            self.assertEqual(len(result), 2)
            for locus in result:
                self.assertFalse(locus.merged_by_ld)
                self.assertEqual(locus.ld_locus_status, "not_evaluated")

    def test_lead_absent_from_vcf_stays_not_evaluated(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            loci = build_loci(tmp, ROWS_A, ROWS_B)
            vcf = tmp / "ld.vcf"
            write_vcf(vcf, [(10, "A", "G", ["0/0", "0/1", "1/1", "1/1"])])
            result = merge_clusters_by_ld(loci, vcf, 0.5)
            self.assertEqual(len(result), 2)
            for locus in result:
                self.assertFalse(locus.merged_by_ld)
                self.assertEqual(locus.ld_locus_status, "not_evaluated")

    def test_transitive_merge_recomputes_lead(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            rows_a = [
                ["s1", "chr1", "10", "3e-8", "A", "G", "0.2"],
                ["s2", "chr1", "80", "1e-8", "C", "T", "0.3"],
                ["s3", "chr1", "150", "2e-8", "G", "A", "0.4"],
            ]
            rows_b = [[r[0], r[1], r[2], "5e-8", r[4], r[5], "0.5"] for r in rows_a]
            loci = build_loci(tmp, rows_a, rows_b)
            self.assertEqual(len(loci), 3)
            vcf = tmp / "ld.vcf"
            gts = ["0/0", "0/1", "1/1", "1/1"]
            write_vcf(vcf, [(10, "A", "G", gts), (80, "C", "T", gts), (150, "G", "A", gts)])
            merged = merge_clusters_by_ld(loci, vcf, 0.5)
            self.assertEqual(len(merged), 1)
            locus = merged[0]
            self.assertEqual(locus.merged_cluster_count, 3)
            self.assertEqual(locus.lead.snp, "s2")  # recomputed min-P lead
            self.assertEqual(locus.start, 10)
            self.assertEqual(locus.end, 150)
            self.assertEqual(len(locus.variants), 3)

    def test_chrom_naming_fallback(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            loci = build_loci(tmp, ROWS_A, ROWS_B)
            vcf = tmp / "ld.vcf"
            write_vcf(
                vcf,
                [
                    (10, "A", "G", ["0/0", "0/1", "1/1", "1/1"]),
                    (80, "C", "T", ["0/0", "0/1", "1/1", "1/1"]),
                ],
                chrom="1",
            )
            merged = merge_clusters_by_ld(loci, vcf, 0.5)
            self.assertEqual(len(merged), 1)
            self.assertTrue(merged[0].merged_by_ld)

    def test_r2_threshold_validation(self):
        loci = (make_locus("chr1", 10, 1e-8),)
        for bad in (0, -0.1, 1.5):
            with self.assertRaisesRegex(ValueError, "r2_threshold"):
                merge_clusters_by_ld(loci, None, bad)


class RankLociTest(unittest.TestCase):
    def test_rank_orders_by_nominal_lead_p(self):
        loci = (
            make_locus("chr1", 300, 1e-3),
            make_locus("chr2", 100, 1e-9),
            make_locus("chr1", 50, 1e-6),
        )
        ranked = rank_loci(loci)
        self.assertEqual([l.lead.pos for l in ranked], [100, 50, 300])
        self.assertEqual([l.locus_rank for l in ranked], [1, 2, 3])
        self.assertEqual([l.priority_peak for l in ranked], [True, False, False])

    def test_rank_ties_break_deterministically(self):
        loci = (make_locus("chr2", 10, 1e-8), make_locus("chr1", 20, 1e-8))
        ranked = rank_loci(loci)
        self.assertEqual([l.chrom for l in ranked], ["chr1", "chr2"])
        self.assertEqual([l.locus_rank for l in ranked], [1, 2])


class PipelineIntegrationTest(unittest.TestCase):
    def test_ld_merge_requires_vcf(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            config = {
                "scenario": "no_vcf", "trait": "height", "data_mode": "full",
                "reference_build": "construct-v1", "ld_merge": True,
                "model_a": {"name": "MLM", "path": "a.csv", "p_threshold": 1e-5},
                "model_b": {"name": "FarmCPU", "path": "b.csv", "p_threshold": 1e-5},
                "gff": "ref.gff", "genome": "genome.fa",
            }
            with self.assertRaisesRegex(ValueError, "ld_merge requires a VCF"):
                run_candidate_pipeline(config, tmp)

    def test_pipeline_merges_ranks_and_reports(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            inputs = tmp / "inputs"
            inputs.mkdir()
            write_signals(inputs / "a.csv", ROWS_A)
            write_signals(inputs / "b.csv", ROWS_B)
            end = len(CODING)
            (inputs / "ref.gff").write_text(
                f"chr1\tconstruct\tgene\t1\t{end}\t.\t+\t.\tID=GeneA\n"
                f"chr1\tconstruct\tmRNA\t1\t{end}\t.\t+\t.\tID=TxA;Parent=GeneA\n"
                f"chr1\tconstruct\tCDS\t1\t{end}\t.\t+\t0\tParent=TxA\n",
                encoding="utf-8",
            )
            (inputs / "genome.fa").write_text(">chr1\n" + CODING + "N" * 37 + "\n", encoding="utf-8")
            write_vcf(inputs / "variants.vcf", [
                (10, "A", "G", ["0/0", "0/1", "1/1", "1/1"]),
                (80, "C", "T", ["0/0", "0/1", "1/1", "1/1"]),
            ])
            config = {
                "scenario": "ld_merge", "base_dir": "inputs", "trait": "height",
                "data_mode": "full", "reference_build": "construct-v1",
                "random_seed": 20261109,
                "model_a": {"name": "MLM", "path": "a.csv", "p_threshold": 1e-5},
                "model_b": {"name": "FarmCPU", "path": "b.csv", "p_threshold": 1e-5},
                "gff": "ref.gff", "genome": "genome.fa", "vcf": "variants.vcf",
                "gap": 20, "flank": 50, "region_flank": 30,
                "ld_merge": True, "ld_merge_r2": 0.5,
            }
            result = run_candidate_pipeline(config, tmp)
            self.assertEqual(len(result["loci"]), 1)
            locus = result["loci"][0]
            self.assertTrue(locus.merged_by_ld)
            self.assertEqual(locus.locus_rank, 1)
            self.assertTrue(locus.priority_peak)
            table = (result["report"]["tables"] / "loci.tsv").read_text(encoding="utf-8")
            self.assertIn("locus_rank", table.splitlines()[0])
            self.assertIn("priority_peak", table.splitlines()[0])
            self.assertIn("merged_by_ld", table.splitlines()[0])
            self.assertEqual(locus.ld_locus_status, "merged_by_ld")
            report_text = result["report"]["markdown"].read_text(encoding="utf-8")
            self.assertIn("nominal lead", report_text)
            self.assertIn("merged_by_ld", report_text)

    def test_pipeline_without_ld_merge_keeps_distance_clusters(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            inputs = tmp / "inputs"
            inputs.mkdir()
            write_signals(inputs / "a.csv", ROWS_A)
            write_signals(inputs / "b.csv", ROWS_B)
            end = len(CODING)
            (inputs / "ref.gff").write_text(
                f"chr1\tconstruct\tgene\t1\t{end}\t.\t+\t.\tID=GeneA\n"
                f"chr1\tconstruct\tmRNA\t1\t{end}\t.\t+\t.\tID=TxA;Parent=GeneA\n"
                f"chr1\tconstruct\tCDS\t1\t{end}\t.\t+\t0\tParent=TxA\n",
                encoding="utf-8",
            )
            (inputs / "genome.fa").write_text(">chr1\n" + CODING + "N" * 37 + "\n", encoding="utf-8")
            write_vcf(inputs / "variants.vcf", [
                (10, "A", "G", ["0/0", "0/1", "1/1", "1/1"]),
                (80, "C", "T", ["0/0", "0/1", "1/1", "1/1"]),
            ])
            config = {
                "scenario": "no_merge", "base_dir": "inputs", "trait": "height",
                "data_mode": "full", "reference_build": "construct-v1",
                "random_seed": 20261109,
                "model_a": {"name": "MLM", "path": "a.csv", "p_threshold": 1e-5},
                "model_b": {"name": "FarmCPU", "path": "b.csv", "p_threshold": 1e-5},
                "gff": "ref.gff", "genome": "genome.fa", "vcf": "variants.vcf",
                "gap": 20, "flank": 50, "region_flank": 30,
            }
            result = run_candidate_pipeline(config, tmp)
            self.assertEqual(len(result["loci"]), 2)
            self.assertTrue(all(not l.merged_by_ld for l in result["loci"]))
            self.assertEqual([l.locus_rank for l in result["loci"]], [1, 2])


if __name__ == "__main__":
    unittest.main()
