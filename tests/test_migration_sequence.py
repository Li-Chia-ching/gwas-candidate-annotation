import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.sequence import (
    ProteinQCError,
    extract_cds_protein,
    extract_proteome,
    parse_gff_models,
    revcomp,
    translate,
    write_sequence_outputs,
)


class SequenceMigrationTest(unittest.TestCase):
    def test_negative_strand_cds_is_reverse_complemented_and_translated(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff = tmp / "ref.gff"
            gff.write_text(
                "chr1\ttest\tgene\t5\t13\t.\t-\t.\tID=GeneNeg\n"
                "chr1\ttest\tmRNA\t5\t13\t.\t-\t.\tID=TxNeg;Parent=GeneNeg\n"
                "chr1\ttest\tCDS\t5\t13\t.\t-\t0\tParent=TxNeg\n",
                encoding="utf-8",
            )
            genome = tmp / "genome.fa"
            genome.write_text(">chr1\nNNNNTTAAGCCATNN\n", encoding="utf-8")
            records = extract_cds_protein(gff, genome, ["GeneNeg"])
            self.assertEqual(records[0].cds_sequence, "ATGGCTTAA")
            self.assertEqual(records[0].protein_sequence, "MA*")
            self.assertEqual(records[0].qc_status, "OK")
            cds_path, protein_path, report_path = write_sequence_outputs(records, tmp / "out")
            self.assertTrue(cds_path.exists())
            self.assertTrue(protein_path.exists())
            self.assertIn("OK", report_path.read_text(encoding="utf-8"))

    def test_bad_protein_qc_blocks_candidate_extraction(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            gff = tmp / "ref.gff"
            gff.write_text(
                "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=BadGene\n"
                "chr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=TxBad;Parent=BadGene\n"
                "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=TxBad\n",
                encoding="utf-8",
            )
            genome = tmp / "genome.fa"
            genome.write_text(">chr1\nATGAAAGGG\n", encoding="utf-8")
            with self.assertRaises(ProteinQCError):
                extract_cds_protein(gff, genome, ["BadGene"])
            proteome = extract_proteome(gff, genome)
            self.assertEqual(proteome[0].qc_status, "FAIL")
            self.assertIn("missing_stop_codon", proteome[0].qc_issues)

    def test_low_level_translation_helpers_are_stable(self):
        self.assertEqual(revcomp("ATGCN"), "NGCAT")
        self.assertEqual(translate("ATGGCTTAA"), "MA*")
        genes, transcripts, cds = parse_gff_models(self._tiny_gff())
        self.assertEqual(genes["GeneA"]["transcripts"], ["TxA"])
        self.assertEqual(transcripts["TxA"], "GeneA")
        self.assertEqual(cds["TxA"], [(1, 9)])

    def _tiny_gff(self):
        tmp = tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8", suffix=".gff")
        with tmp:
            tmp.write(
                "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=GeneA\n"
                "chr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=TxA;Parent=GeneA\n"
                "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=TxA\n"
            )
        self.addCleanup(Path(tmp.name).unlink)
        return Path(tmp.name)


if __name__ == "__main__":
    unittest.main()
