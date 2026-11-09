import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GWAS = ROOT / "quality" / "plant-wgs-gwas-skill" / "scripts" / "03_gwas.sh"


class GemmaInputIntegrationTest(unittest.TestCase):
    def test_gemma_receives_fam_phenotype_numeric_covariates_and_explicit_k(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            work = tmp / "work"
            (work / "qc").mkdir(parents=True)
            (work / "structure").mkdir(parents=True)
            (work / "logs").mkdir(parents=True)
            (work / "qc" / "clean.bed").write_bytes(b"")
            (work / "qc" / "clean.bim").write_text("", encoding="utf-8")
            (work / "qc" / "clean.fam").write_text(
                "S1 S1 0 0 0 -9\nS2 S2 0 0 0 -9\n", encoding="utf-8"
            )
            (work / "structure" / "pca.eigenvec").write_text(
                "S1 S1 0.1 0.2\nS2 S2 -0.1 0.0\n", encoding="utf-8"
            )
            kin = work / "structure" / "kin.cXX.txt"
            kin.write_text("1 0\n0 1\n", encoding="utf-8")
            (work / "structure" / "kinship_gemma.txt").write_text(str(kin), encoding="utf-8")
            phenotype = tmp / "trait.tsv"
            phenotype.write_text("IID\tBLUP_days\nS1\t78.2\nS2\tNA\n", encoding="utf-8")
            fake_gemma = tmp / "fake_gemma.sh"
            fake_gemma.write_text(
                "#!/usr/bin/env bash\n"
                "bfile=''; kfile=''; covar=''; n=''; outdir=''; prefix=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  case \"$1\" in\n"
                "    -bfile) bfile=\"$2\"; shift 2 ;;\n"
                "    -k) kfile=\"$2\"; shift 2 ;;\n"
                "    -c) covar=\"$2\"; shift 2 ;;\n"
                "    -n) n=\"$2\"; shift 2 ;;\n"
                "    -outdir) outdir=\"$2\"; shift 2 ;;\n"
                "    -o) prefix=\"$2\"; shift 2 ;;\n"
                "    *) shift ;;\n"
                "  esac\n"
                "done\n"
                "printf 'bfile=%s\\nk=%s\\nc=%s\\nn=%s\\n' \"$bfile\" \"$kfile\" \"$covar\" \"$n\" > \"$CAPTURE\"\n"
                "mkdir -p \"$outdir\"\n"
                "printf 'chr\\trs\\tps\\tn_miss\\tallele1\\tallele0\\taf\\tbeta\\tse\\tp_wald\\n' > \"$outdir/${prefix}.assoc.txt\"\n",
                encoding="utf-8",
            )
            fake_gemma.chmod(0o755)
            config = tmp / "config.yaml"
            config.write_text(
                "paths:\n"
                "  workdir: ./work\n"
                "genotype:\n"
                "  vcf: genotypes.vcf.gz\n"
                "structure:\n"
                "  n_pc: 2\n"
                "phenotypes:\n"
                "  - name: trait\n"
                "    file: ./trait.tsv\n"
                "    column: BLUP_days\n"
                "    type: quantitative\n"
                "    model: LMM\n"
                "model:\n"
                "  default: LMM\n"
                "software:\n"
                "  plink: plink2\n"
                f"  gemma: {fake_gemma.as_posix()}\n",
                encoding="utf-8",
            )
            capture = tmp / "gemma_args.txt"
            env = dict(os.environ, CAPTURE=str(capture))
            result = subprocess.run(
                ["bash", str(GWAS), str(config)], text=True, capture_output=True, env=env
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            args = dict(
                line.split("=", 1) for line in capture.read_text(encoding="utf-8").splitlines()
            )
            self.assertEqual(args["n"], "1")
            self.assertEqual(Path(args["k"]).name, "kin.cXX.txt")
            self.assertTrue(args["c"].endswith("trait.gemma.covar.txt"))
            self.assertTrue(args["bfile"].endswith("trait_gemma"))
            fam = (work / "assoc" / "trait_gemma.fam").read_text(encoding="utf-8").splitlines()
            self.assertEqual(fam[0].split()[5], "78.2")
            self.assertEqual(fam[1].split()[5], "-9")
            covar_lines = (work / "assoc" / "trait.gemma.covar.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertEqual(covar_lines, ["1\t0.1\t0.2", "1\t-0.1\t0.0"])
            plink_pheno = (work / "assoc" / "trait.pheno").read_text(encoding="utf-8").splitlines()[0]
            self.assertEqual(plink_pheno, "FID\tIID\tBLUP_days")


if __name__ == "__main__":
    unittest.main()
