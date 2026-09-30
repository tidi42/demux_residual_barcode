"""Pin the residual-barcode detector behaviour on constructed reads.

Only NB01 and NB02 (forward + reverse complement) are loaded. F = NB01 forward,
R = NB01 reverse complement, F(1)/R(1) = one substitution at position 10,
numbers = random ACGT insert lengths from a fixed seed.

The advice list has no S2; S2 is defined here as a single terminal copy
(F . 300), the control whose class must not change between detectors.
"""
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import compare_dorado_barbell_outputs as production  # noqa: E402

NB01 = "CACAAAGACACCGACAACTTTCTT"
NB02 = "ACAGACGACTACAAACGGAATCGA"
WINDOW = 80
SEED = 20260929


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def substituted(sequence, position=10):
    replacement = next(base for base in "ACGT" if base != sequence[position])
    return sequence[:position] + replacement + sequence[position + 1:]


F, R = NB01, reverse_complement(NB01)
F1, R1 = substituted(F), substituted(R)
G = NB02

BARCODES = [
    dict(barcode_name=name, barcode_nr=nr, direction=direction, sequence=sequence,
         barcode_length=len(sequence))
    for name, nr, forward in (("NB01", "01", NB01), ("NB02", "02", NB02))
    for direction, sequence in (("forward", forward), ("backward", reverse_complement(forward)))
]

SCENARIOS = {
    "S1": [F, F1, 300, F1],
    "S2": [F, 300],
    "S3": [F, 200, F, 200],
    "S4": [F1, 200, F, 200, F1],
    "S5": [200, F, 200, F],
    "S6": [F, 200, F1, 200],
    "S7": [F1, 200, F, 200],
    "S8": [F1, 200, F1, 200],
    "S9": [F, 200, R, F, 200, R],
    "S10": [F, 200, R1, F1, 200, R],
    "S11": [200, F, 200],
    "S12": [F, 200, G, 200],
}


def build_reads(seed=SEED):
    """Return {name: sequence}; inserts are drawn in scenario order from one generator."""
    generator = random.Random(seed)
    reads = {}
    for name, segments in SCENARIOS.items():
        reads[name] = "".join(
            "".join(generator.choices("ACGT", k=part)) if isinstance(part, int) else part
            for part in segments
        )
    return reads


READS = build_reads()


def retained(scan):
    """(barcode_nr, direction, start) triples parsed from residual_details."""
    if not scan["residual_details"]:
        return []
    triples = []
    for item in scan["residual_details"].split("|"):
        nr, direction, span, _ = item.split(":")
        triples.append((nr, direction, int(span.split("-")[0])))
    return triples


BEST_HIT_EXPECTED = {
    "S1": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S2": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S3": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S4": ("single_residual_barcode__internal", [("01", "forward", 224)]),
    "S5": ("single_residual_barcode__internal", [("01", "forward", 200)]),
    "S6": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S7": ("single_residual_barcode__internal", [("01", "forward", 224)]),
    "S8": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S9": ("multiple_same_barcode_mixed_direction__internal",
           [("01", "forward", 0), ("01", "backward", 224)]),
    "S10": ("multiple_same_barcode_mixed_direction__both_ends",
            [("01", "forward", 0), ("01", "backward", 472)]),
    "S11": ("single_residual_barcode__internal", [("01", "forward", 200)]),
    "S12": ("multiple_different_barcodes__internal", [("01", "forward", 0), ("02", "forward", 224)]),
}


class BestHitBehaviourTests(unittest.TestCase):
    """Step 0: the legacy one-best-match-per-pattern detector, pinned before any change."""

    def setUp(self):
        production.init_worker(BARCODES, 1, 1, 1, "best_hit")

    def test_scenarios(self):
        for name, (pattern_class, triples) in BEST_HIT_EXPECTED.items():
            with self.subTest(read=name):
                scan = production.scan_record(READS[name], "NB01", WINDOW)
                self.assertEqual(scan["pattern_class"], pattern_class)
                self.assertEqual(retained(scan), triples)


SAME = "multiple_same_barcode_same_direction"
MIXED = "multiple_same_barcode_mixed_direction"
OCCURRENCES_EXPECTED = {
    "S1": (SAME + "__both_ends", [("01", "forward", 0), ("01", "forward", 24), ("01", "forward", 348)]),
    "S2": ("single_residual_barcode__terminal", [("01", "forward", 0)]),
    "S3": (SAME + "__internal", [("01", "forward", 0), ("01", "forward", 224)]),
    "S4": (SAME + "__internal", [("01", "forward", 0), ("01", "forward", 224), ("01", "forward", 448)]),
    "S5": (SAME + "__internal", [("01", "forward", 200), ("01", "forward", 424)]),
    "S6": (SAME + "__internal", [("01", "forward", 0), ("01", "forward", 224)]),
    "S7": (SAME + "__internal", [("01", "forward", 0), ("01", "forward", 224)]),
    "S8": (SAME + "__internal", [("01", "forward", 0), ("01", "forward", 224)]),
    "S9": (MIXED + "__internal", [("01", "forward", 0), ("01", "backward", 224),
                                  ("01", "forward", 248), ("01", "backward", 472)]),
    "S10": (MIXED + "__internal", [("01", "forward", 0), ("01", "backward", 224),
                                   ("01", "forward", 248), ("01", "backward", 472)]),
    "S11": ("single_residual_barcode__internal", [("01", "forward", 200)]),
    "S12": ("multiple_different_barcodes__internal", [("01", "forward", 0), ("02", "forward", 224)]),
}
# Reads whose binary internal status differs between the detectors.
NEWLY_INTERNAL = {"S3", "S6", "S8", "S10"}


class OccurrenceDetectorTests(unittest.TestCase):
    """Step 4: the same reads through --detector occurrences."""

    def setUp(self):
        production.init_worker(BARCODES, 1, 1, 1, "occurrences")

    def tearDown(self):
        production.init_worker(BARCODES, 1, 1, 1, "best_hit")

    def test_scenarios(self):
        for name, (pattern_class, triples) in OCCURRENCES_EXPECTED.items():
            with self.subTest(read=name):
                scan = production.scan_record(READS[name], "NB01", WINDOW)
                self.assertEqual(scan["pattern_class"], pattern_class)
                self.assertEqual(retained(scan), triples)
                self.assertEqual(scan["residual_count"], len(triples))

    def test_s1_counts_three_occurrences(self):
        self.assertEqual(production.scan_record(READS["S1"], "NB01", WINDOW)["residual_count"], 3)

    def test_internal_status_changes_only_for_masked_reads(self):
        changed = {name for name in SCENARIOS
                   if BEST_HIT_EXPECTED[name][0].endswith("__internal")
                   != OCCURRENCES_EXPECTED[name][0].endswith("__internal")}
        self.assertEqual(changed, NEWLY_INTERNAL)

    def test_barcode_names_are_mapped(self):
        production.init_worker(BARCODES, 1, 1, 1, "occurrences")
        matches = production.occurrence_matches(READS["S12"], production.WORKER_SCANNER,
                                                production.WORKER_BARCODE_NAMES)
        self.assertEqual([m["barcode_name"] for m in matches], ["NB01", "NB02"])
        self.assertEqual(matches[1]["matched_sequence"], NB02)


def mutate(sequence, generator):
    """Up to one substitution, one insertion and one deletion (the production edit budget)."""
    bases = list(sequence)
    if generator.random() < 0.7:
        i = generator.randrange(len(bases))
        bases[i] = generator.choice([b for b in "ACGT" if b != bases[i]])
    if generator.random() < 0.5:
        bases.insert(generator.randrange(len(bases) + 1), generator.choice("ACGT"))
    if generator.random() < 0.5:
        del bases[generator.randrange(len(bases))]
    return "".join(bases)


class DetectorAgreementTests(unittest.TestCase):
    """A read with no best-hit match has no occurrence match, and vice versa."""

    def test_random_reads_agree_on_no_match(self):
        generator = random.Random(SEED + 1)
        reads = []
        for index in range(1000):
            read = "".join(generator.choices("ACGT", k=generator.randrange(30, 600)))
            if index % 4 == 0:
                barcode = mutate(generator.choice(BARCODES)["sequence"], generator)
                position = generator.randrange(len(read) + 1)
                read = read[:position] + barcode + read[position:]
            reads.append(read)
        positive = {}
        for detector in ("best_hit", "occurrences"):
            production.init_worker(BARCODES, 1, 1, 1, detector)
            positive[detector] = [production.scan_record(read, "NB01", WINDOW)["residual_count"] > 0
                                  for read in reads]
        production.init_worker(BARCODES, 1, 1, 1, "best_hit")
        self.assertEqual(positive["best_hit"], positive["occurrences"])
        self.assertGreater(sum(positive["best_hit"]), 200)
        self.assertLess(sum(positive["best_hit"]), 1000)


class DecoyTests(unittest.TestCase):
    def test_decoys_preserve_composition_and_pairing(self):
        sets = production.make_decoy_sets(BARCODES, 2, 7)
        self.assertEqual(len(sets), 2)
        self.assertEqual(sets, production.make_decoy_sets(BARCODES, 2, 7))
        for decoys in sets:
            by_key = {(d["barcode_nr"], d["direction"]): d["sequence"] for d in decoys}
            for real in BARCODES:
                decoy = by_key[(real["barcode_nr"], real["direction"])]
                self.assertEqual(sorted(decoy), sorted(real["sequence"]))
                self.assertNotEqual(decoy, real["sequence"])
            for nr in ("01", "02"):
                self.assertEqual(by_key[(nr, "backward")], reverse_complement(by_key[(nr, "forward")]))


if __name__ == "__main__":
    unittest.main()
