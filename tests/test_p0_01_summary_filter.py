import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "legacy" / "gwas-candidate-annotation" / "scripts" / "screen_dual_model_leads.py"


class SummarySignificanceFilterTest(unittest.TestCase):
    def test_full_summary_is_filtered_before_model_intersection(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff = tmp / "ref.gff"
            gff.write_text("chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=GeneA\n", encoding="utf-8")
            a = tmp / "a.csv"
            b = tmp / "b.csv"
            header = ["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect"]
            a_rows = [
                ["s1", "chr1", "100", "1e-6", "A", "G", "0.2"],
                ["s2", "chr1", "200", "1e-9", "C", "T", "0.4"],
                ["s3", "chr1", "300", "1e-3", "G", "A", "0.1"],
                ["a_only", "chr1", "400", "1e-9", "T", "C", "0.2"],
            ]
            b_rows = [
                ["s1", "chr1", "100", "2e-6", "A", "G", "0.3"],
                ["s2", "chr1", "200", "5e-8", "C", "T", "0.5"],
                ["s3", "chr1", "300", "4e-8", "G", "A", "0.2"],
                ["b_only", "chr1", "500", "1e-9", "T", "A", "0.1"],
            ]
            for path, rows in ((a, a_rows), (b, b_rows)):
                with path.open("w", newline="", encoding="utf-8") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(header)
                    writer.writerows(rows)
            outdir = tmp / "out"
            cmd = [
                sys.executable, str(SCRIPT),
                "--signals-a", str(a), "--signals-b", str(b),
                "--model-a", "MLM", "--model-b", "FarmCPU",
                "--gff", str(gff), "--trait", "T", "--outdir", str(outdir),
                "--p-threshold-a", "1e-5", "--p-threshold-b", "1e-5",
                "--reference-build", "buildA", "--gff-reference-build", "buildA",
                "--data-mode", "full",
            ]
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with (outdir / "T_dual_model_snps.csv").open(newline="", encoding="utf-8-sig") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual([row["pos"] for row in rows], ["100", "200"])


if __name__ == "__main__":
    unittest.main()
