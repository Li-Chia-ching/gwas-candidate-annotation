import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "legacy" / "gwas-candidate-annotation" / "scripts" / "screen_dual_model_leads.py"


class DistanceClusterSeparationTest(unittest.TestCase):
    def test_distance_cluster_is_not_reported_as_ld_locus(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff = tmp / "ref.gff"
            gff.write_text("chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=GeneA\n", encoding="utf-8")
            a = tmp / "a.csv"
            b = tmp / "b.csv"
            for path in (a, b):
                with path.open("w", newline="", encoding="utf-8") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(["SNP", "CHROM", "POS", "P.value", "REF", "ALT", "Effect"])
                    writer.writerows([
                        ["s1", "chr1", "100", "1e-8", "A", "G", "0.2"],
                        ["s2", "chr1", "150", "2e-8", "C", "T", "0.3"],
                    ])
            outdir = tmp / "out"
            cmd = [
                sys.executable, str(SCRIPT),
                "--signals-a", str(a), "--signals-b", str(b),
                "--model-a", "MLM", "--model-b", "FarmCPU",
                "--gff", str(gff), "--trait", "T", "--outdir", str(outdir),
                "--p-threshold-a", "1e-5", "--p-threshold-b", "1e-5",
                "--reference-build", "buildA", "--gff-reference-build", "buildA",
                "--data-mode", "full",
                "--gap", "100",
            ]
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertNotIn("independent loci", result.stdout.lower())
            with (outdir / "T_dual_model_leads.csv").open(newline="", encoding="utf-8-sig") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["cluster_type"], "distance_cluster")
            self.assertEqual(rows[0]["ld_locus_status"], "not_evaluated")


if __name__ == "__main__":
    unittest.main()
