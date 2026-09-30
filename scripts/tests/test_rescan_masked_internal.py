"""rescan_masked_internal.py on a small database built with production best_hit scan_record()."""
import csv
import sqlite3
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import compare_dorado_barbell_outputs as production  # noqa: E402
import rescan_masked_internal  # noqa: E402
from test_best_hit_behaviour import BARCODES, NB01, NB02, NEWLY_INTERNAL, READS, WINDOW  # noqa: E402


class RescanMaskedInternalTests(unittest.TestCase):
    def test_reclassifies_exactly_masked_reads(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "old_host" / "c4" / "NB01.fastq"
            fastq.parent.mkdir(parents=True)
            with fastq.open("w") as handle:
                for name, sequence in READS.items():
                    handle.write(f"@{name} extra\n{sequence}\n+\n{'I' * len(sequence)}\n")
                handle.write(f"@empty\n{'ACGT' * 50}\n+\n{'I' * 200}\n")
            csv_path = root / "barcodes.csv"
            csv_path.write_text("Component,Forward sequence,Reverse sequence\n"
                                f"NB01,{NB01},{production.reverse_complement(NB01)}\n"
                                f"NB02,{NB02},{production.reverse_complement(NB02)}\n")
            self.assertEqual(production.load_barcodes(csv_path), BARCODES)

            db = root / "cmp_all" / "comparison.sqlite"
            db.parent.mkdir()
            con = production.connect_db(db)
            production.init_db(con)
            production.init_worker(BARCODES, 1, 1, 1, "best_hit")
            rows, _ = production.process_fastq_range(("c4", str(fastq), "01", WINDOW, 0, len(READS) + 1))
            production.insert_read_rows(con, rows)
            con.commit()
            con.close()

            # The database records the old path; the FASTQ now lives under new_host.
            moved = root / "new_host"
            (root / "old_host").rename(moved)
            out = root / "rescan.csv"
            reclassified = root / "reclassified_LOCAL.csv"
            code = rescan_masked_internal.main([
                "--db", str(db), "--barcode-csv", str(csv_path), "--out", str(out),
                "--threads", "1", "--terminal-window", str(WINDOW), "--experiment", "test",
                "--path-map", f"{root / 'old_host'}={moved}", "--reclassified-out", str(reclassified),
            ])
            self.assertEqual(code, 0)
            with sqlite3.connect(db) as check:
                self.assertEqual(check.execute("SELECT COUNT(*) FROM reads").fetchone()[0], len(READS) + 1)

            [row] = list(csv.DictReader(out.read_text().splitlines()))
            changed = {(r["read_id"], r["change"]) for r in csv.DictReader(reclassified.read_text().splitlines())}
            self.assertEqual(changed, {(name, "newly_internal") for name in NEWLY_INTERNAL})
            legacy_internal = sum(1 for name in READS if name in {"S4", "S5", "S7", "S9", "S11", "S12"})
            self.assertEqual(int(row["total_reads"]), len(READS) + 1)
            self.assertEqual(int(row["legacy_internal_reads"]), legacy_internal)
            self.assertEqual(int(row["rescanned_reads"]), len(READS) - legacy_internal)
            self.assertEqual(int(row["reads_found"]), len(READS) - legacy_internal)
            self.assertEqual(int(row["newly_internal_reads"]), len(NEWLY_INTERNAL))
            self.assertEqual(int(row["internal_reads_lost"]), 0)
            self.assertEqual(int(row["corrected_internal_reads"]), legacy_internal + len(NEWLY_INTERNAL))
            self.assertEqual(int(row["missing_reads"]), 0)
            self.assertAlmostEqual(float(row["corrected_internal_RPM"]),
                                   1e6 * (legacy_internal + len(NEWLY_INTERNAL)) / (len(READS) + 1))

    def test_split_tasks_partitions_reads(self):
        selection = {"a.fastq": {f"r{i}": [("c0", 10, False)] for i in range(10)}}
        tasks = rescan_masked_internal.split_tasks(selection, 3)
        self.assertEqual(len(tasks), 4)
        merged = {}
        for path, share in tasks:
            self.assertEqual(path, "a.fastq")
            self.assertFalse(set(share) & set(merged))
            merged.update(share)
        self.assertEqual(merged, selection["a.fastq"])


if __name__ == "__main__":
    unittest.main()
