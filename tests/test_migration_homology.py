import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.homology import (
    parse_blast_hits,
    read_db_config,
    self_blast_candidates,
    top_homology_candidates,
    write_homology_tables,
)


class HomologyMigrationTest(unittest.TestCase):
    def test_blast_hits_are_homology_candidates_not_orthologs(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            blast = tmp / "Self.blast.tsv"
            blast.write_text(
                "GeneA|Tx1\tGeneA\t100.0\t90\t0\t0\t1\t90\t1\t90\t1e-50\t180\t100\t90\t90\n"
                "GeneA|Tx1\tGeneB\t96.0\t88\t2\t0\t1\t88\t1\t88\t1e-40\t160\t97.8\t90\t90\n",
                encoding="utf-8",
            )
            hits = parse_blast_hits(blast, "Self")
            self.assertEqual(hits[0].evidence_label, "homology_candidate")
            self.assertEqual(hits[0].subject_cover_pct, 100.0)
            top = top_homology_candidates(hits, top=2)
            self.assertEqual(len(top), 2)
            self_hits = self_blast_candidates(top, "Self")
            self.assertEqual(self_hits[1].evidence_label, "homoeolog_candidate")
            summary, self_path = write_homology_tables(hits, tmp / "out", self_database="Self")
            text = summary.read_text(encoding="utf-8") + self_path.read_text(encoding="utf-8")
            self.assertIn("homology_candidate", text)
            self.assertIn("homoeolog_candidate", text)
            self.assertNotIn("ortholog", text.lower())
            self.assertNotIn("copy_number", text.lower())

    def test_db_config_and_zero_hit_boundaries(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            db = tmp / "dbs.tsv"
            db.write_text("Self\tself.pep.fa\n# comment\n", encoding="utf-8")
            self.assertEqual(read_db_config(db), [("Self", tmp / "self.pep.fa")])
            empty = tmp / "empty.blast.tsv"
            empty.write_text("", encoding="utf-8")
            self.assertEqual(parse_blast_hits(empty, "Self"), [])
            self.assertEqual(self_blast_candidates([], "Self"), [])


if __name__ == "__main__":
    unittest.main()
