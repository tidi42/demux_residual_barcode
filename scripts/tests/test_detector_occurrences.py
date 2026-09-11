import importlib.util
import random
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validate_detector_occurrences import OccurrenceScanner, containing_record, internal_count, legacy_matches


class OccurrenceTests(unittest.TestCase):
    def setUp(self):
        self.barcode = "ACGTGACCTGATCGTAGCTACGTA"
        self.scanner = OccurrenceScanner([
            dict(barcode_nr="01", direction="forward", sequence=self.barcode),
        ])

    def test_terminal_and_internal_copies(self):
        sequence = self.barcode + "N" * 200 + self.barcode + "N" * 200
        matches = self.scanner.scan(sequence)
        self.assertEqual([(record["start"], record["end"]) for record in matches], [(0, 24), (224, 248)])
        self.assertEqual(internal_count(matches, len(sequence), 80), 1)
        self.assertEqual(self.scanner.patterns[0].search(sequence).start(), 0)

    def test_approximate_copy_before_exact_copy(self):
        mutated = self.barcode[:8] + "A" + self.barcode[9:]
        sequence = "N" * 10 + mutated + "N" * 200 + self.barcode + "N" * 200
        matches = self.scanner.scan(sequence)
        self.assertEqual(len(matches), 2)
        self.assertEqual([record["edit_distance"] for record in matches], [1, 0])

    def test_overlap_policy_prefers_exact_and_keeps_adjacent(self):
        matches = self.scanner.scan(self.barcode + self.barcode)
        self.assertEqual([(record["start"], record["end"]) for record in matches], [(0, 24), (24, 48)])
        self.assertEqual(self.scanner.scan("N" * 80), [])

    def test_window_boundary_and_short_read(self):
        match = [dict(start=81, end=105)]
        self.assertEqual(internal_count(match, 186, 80), 1)
        self.assertEqual(internal_count(match, 185, 80), 0)
        self.assertEqual(internal_count([dict(start=20, end=44)], 64, 80), 0)

    def test_seed_candidates_equal_brute_force(self):
        generator = random.Random(20260911)
        for trial in range(20):
            altered = list(self.barcode)
            substitution = generator.randrange(24)
            altered[substitution] = generator.choice([base for base in "ACGT" if base != altered[substitution]])
            altered.insert(generator.randrange(25), generator.choice("ACGT"))
            del altered[generator.randrange(len(altered))]
            sequence = "".join(generator.choices("ACGT", k=15)) + "".join(altered) + "N" * 15
            expected = set()
            for start in range(len(sequence)):
                for length in (23, 24, 25):
                    if start + length <= len(sequence) and self.scanner.patterns[0].fullmatch(sequence[start:start + length]):
                        expected.add((start, start + length))
            actual = {(record["start"], record["end"]) for record in self.scanner.candidates(sequence)}
            with self.subTest(trial=trial):
                self.assertEqual(actual, expected)

    def test_random_byte_resolves_containing_record(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            records = [b"@record1\nACGT\n+\n@@@@\n", b"@record2\nACGTAC\n+\nIIIIII\n"]
            path.write_bytes(b"".join(records))
            for position in range(sum(map(len, records))):
                expected = (0, len(records[0]), "ACGT") if position < len(records[0]) else (
                    len(records[0]), len(records[1]), "ACGTAC")
                self.assertEqual(containing_record(path, position), expected)

    def test_legacy_arm_matches_production_scanner(self):
        path = Path(__file__).resolve().parents[1] / "compare_dorado_barbell_outputs.py"
        specification = importlib.util.spec_from_file_location("production_scanner", path)
        production = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(production)
        barcodes = [dict(barcode_name="NB01", barcode_nr="01", direction="forward", sequence=self.barcode)]
        production.init_worker(barcodes, 1, 1, 1)
        for sequence in [self.barcode + "N"*200 + self.barcode + "N"*200,
                         "N"*100 + self.barcode + "N"*100, "N"*100]:
            for window in (40, 80, 120):
                matches = legacy_matches(self.scanner, sequence)
                result = production.scan_record(sequence, "NB01", window)
                self.assertEqual(result["residual_count"], len(matches))
                self.assertEqual(result["internal_residual_count"], internal_count(matches, len(sequence), window))


if __name__ == "__main__":
    unittest.main()