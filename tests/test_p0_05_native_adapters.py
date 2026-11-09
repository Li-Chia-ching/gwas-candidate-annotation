import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "04_plots.py"
CHECKS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "05_checks.py"


class NativeAssociationAdapterTest(unittest.TestCase):
    def run_tools(self, path, tmp):
        plot_out = tmp / "figures"
        report = tmp / "report.md"
        plot = subprocess.run([
            sys.executable, str(PLOTS), "--input", str(path), "--pheno", "T",
            "--out", str(plot_out), "--threshold", "1e-3", "--data-mode", "full",
        ], text=True, capture_output=True)
        check = subprocess.run([
            sys.executable, str(CHECKS), "--input", str(path), "--pheno", "T",
            "--out", str(report), "--threshold", "1e-3", "--expected-n", "250", "--data-mode", "full",
        ], text=True, capture_output=True)
        return plot, check, plot_out, report

    def test_plink2_native_columns_are_adapted(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            path = tmp / "plink2.tsv"
            path.write_text(
                "#CHROM\tPOS\tID\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tP\n"
                "1\t100\trs1\tA\tG\tG\t250\t0.2\t0.05\t1e-6\n"
                "1\t200\trs2\tC\tT\tT\t250\t0.1\t0.05\t0.02\n",
                encoding="utf-8",
            )
            plot, check, plot_out, report = self.run_tools(path, tmp)
            self.assertEqual(plot.returncode, 0, plot.stderr + plot.stdout)
            self.assertEqual(check.returncode, 0, check.stderr + check.stdout)
            self.assertTrue((plot_out / "manhattan.png").exists())
            self.assertTrue((plot_out / "qq.png").exists())
            self.assertTrue(report.exists())

    def test_gemma_native_columns_are_adapted(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            path = tmp / "gemma.tsv"
            path.write_text(
                "chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tp_wald\n"
                "1\trs1\t100\t0\tG\tA\t0.4\t0.2\t0.05\t1e-6\n"
                "1\trs2\t200\t0\tT\tC\t0.3\t0.1\t0.05\t0.02\n",
                encoding="utf-8",
            )
            plot, check, plot_out, report = self.run_tools(path, tmp)
            self.assertEqual(plot.returncode, 0, plot.stderr + plot.stdout)
            self.assertEqual(check.returncode, 0, check.stderr + check.stdout)
            self.assertTrue((plot_out / "manhattan.png").exists())
            self.assertTrue((plot_out / "qq.png").exists())
            self.assertTrue(report.exists())


if __name__ == "__main__":
    unittest.main()
