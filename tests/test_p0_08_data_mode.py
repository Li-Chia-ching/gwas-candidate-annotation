import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "04_plots.py"
CHECKS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "05_checks.py"
SCREEN = ROOT / "legacy" / "gwas-candidate-annotation" / "scripts" / "screen_dual_model_leads.py"


class DataModeGateTest(unittest.TestCase):
    def make_assoc(self, tmp):
        path = tmp / "assoc.tsv"
        path.write_text(
            "CHR\tBP\tP\tBETA\tN\n"
            "1\t100\t1e-8\t0.2\t250\n"
            "1\t200\t0.2\t0.1\t250\n",
            encoding="utf-8",
        )
        return path

    def run_tools(self, tmp, mode):
        assoc = self.make_assoc(tmp)
        plot = subprocess.run([
            sys.executable, str(PLOTS), "--input", str(assoc), "--pheno", "T",
            "--out", str(tmp / "figures"), "--threshold", "1e-6", "--data-mode", mode,
        ], text=True, capture_output=True)
        check = subprocess.run([
            sys.executable, str(CHECKS), "--input", str(assoc), "--pheno", "T",
            "--out", str(tmp / "report.md"), "--threshold", "1e-6", "--data-mode", mode,
        ], text=True, capture_output=True)
        return plot, check

    def test_full_mode_generates_qq_and_lambda(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            plot, check = self.run_tools(tmp, "full")
            self.assertEqual(plot.returncode, 0, plot.stderr + plot.stdout)
            self.assertEqual(check.returncode, 0, check.stderr + check.stdout)
            self.assertTrue((tmp / "figures" / "qq.png").exists())
            self.assertIn("data_mode=full", plot.stdout)
            self.assertIn("data_mode=full", (tmp / "report.md").read_text(encoding="utf-8"))

    def test_truncated_mode_disables_qq_and_lambda(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            plot, check = self.run_tools(tmp, "truncated")
            self.assertEqual(plot.returncode, 0, plot.stderr + plot.stdout)
            self.assertEqual(check.returncode, 0, check.stderr + check.stdout)
            self.assertFalse((tmp / "figures" / "qq.png").exists())
            self.assertIn("lambda=NA", plot.stdout)
            report = (tmp / "report.md").read_text(encoding="utf-8")
            self.assertIn("data_mode=truncated", report)
            self.assertIn("not evaluated for truncated", report)

    def test_candidate_screen_records_truncated_mode(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            (tmp / "ref.gff").write_text(
                "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=GeneA\n", encoding="utf-8"
            )
            paths = []
            for name in ("a", "b"):
                path = tmp / f"{name}.csv"
                path.write_text(
                    "SNP,CHROM,POS,P.value,REF,ALT,Effect\n"
                    "s1,chr1,100,1e-8,A,G,0.2\n",
                    encoding="utf-8",
                )
                paths.append(path)
            result = subprocess.run([
                sys.executable, str(SCREEN),
                "--signals-a", str(paths[0]), "--signals-b", str(paths[1]),
                "--model-a", "MLM", "--model-b", "FarmCPU",
                "--gff", str(tmp / "ref.gff"), "--trait", "T", "--outdir", str(tmp / "out"),
                "--p-threshold-a", "1e-5", "--p-threshold-b", "1e-5",
                "--reference-build", "buildA", "--gff-reference-build", "buildA",
                "--data-mode", "truncated",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with (tmp / "out" / "T_dual_model_leads.csv").open(newline="", encoding="utf-8-sig") as fh:
                row = list(csv.DictReader(fh))[0]
            self.assertEqual(row["data_mode"], "truncated")


if __name__ == "__main__":
    unittest.main()
