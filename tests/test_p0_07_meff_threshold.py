import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "04_plots.py"
CHECKS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "05_checks.py"


class MeffThresholdConsistencyTest(unittest.TestCase):
    def test_plots_and_checks_use_same_bonferroni_meff_threshold(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            work = tmp / "work"
            (work / "structure").mkdir(parents=True)
            (work / "structure" / "pruned.prune.in").write_text(
                "rs1\nrs2\nrs3\nrs4\n", encoding="utf-8"
            )
            assoc = tmp / "assoc.tsv"
            assoc.write_text(
                "CHR\tBP\tP\tBETA\tN\n"
                "1\t100\t0.001\t0.2\t250\n"
                "1\t200\t0.020\t0.1\t250\n",
                encoding="utf-8",
            )
            config = tmp / "config.yaml"
            config.write_text(
                "paths:\n"
                "  workdir: ./work\n"
                "population:\n"
                "  n_samples_expected: 250\n"
                "threshold:\n"
                "  method: bonferroni_meff\n"
                "  value: 5.0e-8\n"
                "  meff_source: ld_pruned_marker_count\n"
                "diagnostics:\n"
                "  lambda_warn: 1.2\n",
                encoding="utf-8",
            )
            plot = subprocess.run([
                sys.executable, str(PLOTS), "--input", str(assoc), "--pheno", "T",
                "--out", str(tmp / "figures"), "--config", str(config), "--data-mode", "full",
            ], text=True, capture_output=True)
            self.assertEqual(plot.returncode, 0, plot.stderr + plot.stdout)
            check = subprocess.run([
                sys.executable, str(CHECKS), "--input", str(assoc), "--pheno", "T",
                "--out", str(tmp / "report.md"), "--config", str(config), "--data-mode", "full",
            ], text=True, capture_output=True)
            self.assertEqual(check.returncode, 0, check.stderr + check.stdout)
            expected = 0.05 / 4
            self.assertIn(f"threshold={expected:.6g}", plot.stdout)
            report = (tmp / "report.md").read_text(encoding="utf-8")
            self.assertIn(f"threshold={expected:.6g}", report)
            self.assertIn("bonferroni_meff", plot.stdout)
            self.assertIn("bonferroni_meff", report)


if __name__ == "__main__":
    unittest.main()
