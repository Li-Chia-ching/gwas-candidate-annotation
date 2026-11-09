import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "legacy" / "gwas-candidate-annotation" / "scripts" / "screen_dual_model_leads.py"


class AlleleEffectReferenceGateTest(unittest.TestCase):
    def make_inputs(self, tmp, b_ref="A", b_alt="G", effect_a="0.2", effect_b="0.3"):
        gff = tmp / "ref.gff"
        gff.write_text("chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=GeneA\n", encoding="utf-8")
        a = tmp / "a.csv"
        b = tmp / "b.csv"
        with a.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect"])
            writer.writerow(["s1", "chr1", "100", "1e-8", "A", "G", effect_a])
        with b.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect"])
            writer.writerow(["s1", "chr1", "100", "2e-8", b_ref, b_alt, effect_b])
        return gff, a, b

    def run_screen(self, tmp, gff, a, b, reference_build="buildA", gff_build="buildA"):
        return subprocess.run([
            sys.executable, str(SCRIPT),
            "--signals-a", str(a), "--signals-b", str(b),
            "--model-a", "MLM", "--model-b", "FarmCPU",
            "--gff", str(gff), "--trait", "T", "--outdir", str(tmp / "out"),
            "--p-threshold-a", "1e-5", "--p-threshold-b", "1e-5",
            "--reference-build", reference_build,
            "--gff-reference-build", gff_build,
            "--data-mode", "full",
        ], text=True, capture_output=True)

    def test_matching_alleles_effects_and_reference_pass(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff, a, b = self.make_inputs(tmp)
            result = self.run_screen(tmp, gff, a, b)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with (tmp / "out" / "T_dual_model_leads.csv").open(newline="", encoding="utf-8-sig") as fh:
                row = list(csv.DictReader(fh))[0]
            self.assertEqual(row["ref"], "A")
            self.assertEqual(row["alt"], "G")
            self.assertEqual(row["reference_build"], "buildA")
            self.assertEqual(row["effect_direction"], "consistent")

    def test_conflicting_alleles_block(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff, a, b = self.make_inputs(tmp, b_ref="A", b_alt="T")
            result = self.run_screen(tmp, gff, a, b)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("allele mismatch", result.stderr + result.stdout)

    def test_opposite_effect_direction_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff, a, b = self.make_inputs(tmp, effect_b="-0.3")
            result = self.run_screen(tmp, gff, a, b)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("opposite effect directions", result.stderr + result.stdout)

    def test_reference_build_mismatch_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff, a, b = self.make_inputs(tmp)
            result = self.run_screen(tmp, gff, a, b, gff_build="buildB")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("reference build mismatch", result.stderr + result.stdout)


if __name__ == "__main__":
    unittest.main()
