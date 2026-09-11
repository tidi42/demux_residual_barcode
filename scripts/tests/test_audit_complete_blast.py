import sys
import unittest
from tempfile import NamedTemporaryFile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_complete_blast import (
    complete_identity, element_category, historical_identity, read_hsps, taxonomy_key,
)


class CompleteIdentityTests(unittest.TestCase):
    def setUp(self):
        self.record = dict(qlen=24, pident=100, mismatch=0, gapopen=0,
                           qcovs=100, length=24, qstart=1, qend=24,
                           sstart=101, send=124)

    def test_complete_both_coordinate_orientations(self):
        for query_start, query_end in [(1, 24), (24, 1)]:
            for subject_start, subject_end in [(101, 124), (124, 101)]:
                with self.subTest(query_start=query_start, subject_start=subject_start):
                    record = dict(self.record, qstart=query_start, qend=query_end,
                                  sstart=subject_start, send=subject_end)
                    self.assertTrue(complete_identity(record))

    def test_subject_coverage_does_not_prove_complete_hsp(self):
        record = dict(self.record, qstart=2, length=23, send=123)
        self.assertTrue(historical_identity(record))
        self.assertFalse(complete_identity(record))

    def test_each_constraint_is_required(self):
        for change in [dict(pident=99), dict(mismatch=1), dict(gapopen=1),
                       dict(gaps=1), dict(qlen=0), dict(qstart=2), dict(qend=23),
                       dict(length=25), dict(send=125), dict(sstart=0, send=23)]:
            with self.subTest(change=change):
                self.assertFalse(complete_identity(dict(self.record, **change)))

    def test_adapter_query_length(self):
        self.assertTrue(complete_identity(dict(self.record, qlen=47, length=47,
                                              qend=47, send=147)))
        self.assertFalse(complete_identity(dict(self.record, qlen=47)))

    def test_query_coverage_is_not_the_controlling_field(self):
        self.assertTrue(complete_identity(dict(self.record, qcovs=99)))

    def test_exported_reverse_identifier(self):
        with NamedTemporaryFile(mode="w+", suffix=".tsv") as stream:
            stream.write("NB01_reverse\tACC.1\t123\tN/A\t100\t24\t24\t100\t0\t0"
                         "\t1\t24\t124\t101\t0.01\t48.1\tRecord title\n")
            stream.flush()
            record = next(read_hsps(stream.name))
            self.assertEqual(record["orientation"], "revcomp")
            self.assertTrue(complete_identity(record))

    def test_taxonomy_does_not_parse_titles_or_split_ambiguous_sets(self):
        self.assertEqual(taxonomy_key("123"), "taxid:123")
        self.assertEqual(taxonomy_key("456;123;123"), "taxid_set:123;456")
        self.assertEqual(taxonomy_key("N/A"), "unresolved")

    def test_category_precedence_and_neutral_genome_label(self):
        self.assertEqual(element_category("Phage plasmid complete genome"), "Plasmid")
        self.assertEqual(element_category("Eukaryotic chromosome"), "Chromosome/genome")
        self.assertEqual(element_category("Unknown record"), "Other/uncultured")


if __name__ == "__main__":
    unittest.main()