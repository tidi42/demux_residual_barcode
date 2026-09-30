#!/usr/bin/env python3
"""Correct the internal-residual endpoint of a finished best_hit run without a full rescan.

The historical best_hit detector keeps one match per barcode orientation (fewest edits,
ties to the most 5' copy), so an internal copy can be hidden behind a terminal copy of the
same barcode orientation. Only reads with at least one best-hit match can hide a copy: a read
without any best-hit match has no occurrence either. This script therefore rescans only the
reads whose pattern_class is neither no_residual_barcode, no_residual_scan nor *__internal
with OccurrenceScanner, and counts a read as newly internal if
internal_count(matches, len(seq), window) > 0.

Reads already classed *__internal are rechecked as well (they are few). Greedy overlap
resolution can, in rare cases, replace an internal best-hit match with an overlapping
occurrence that crosses the terminal-window boundary; such reads are reported as
internal_reads_lost. With that recheck, corrected_internal_reads equals the internal read
count of a full --detector occurrences rescan of the same FASTQs.

Output: one CSV row per tool (condition). No read identifiers are written to the CSV;
--reclassified-out optionally lists them for local inspection (keep that file private).

Example:
  python3 scripts/rescan_masked_internal.py \\
    --db outputs/experiment_a_demux/cmp_all/comparison.sqlite \\
    --barcode-csv barcodes.csv --experiment experiment_a \\
    --path-map /old/data/root=/new/data/root \\
    --out rescan/experiment_a_rescan_masked_internal.csv --threads 16
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
import time
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_dorado_barbell_outputs import extract_read_id, fastq_reader, load_barcodes  # noqa: E402
from validate_detector_occurrences import OccurrenceScanner, internal_count  # noqa: E402

NOT_RESCANNED = ("no_residual_barcode", "no_residual_scan")
FIELDS = [
    "experiment", "tool", "total_reads", "legacy_internal_reads", "legacy_internal_RPM",
    "rescanned_reads", "reads_found", "newly_internal_reads", "fraction_of_rescanned_reclassified",
    "internal_reads_rechecked", "internal_reads_found", "internal_reads_lost",
    "corrected_internal_reads", "corrected_internal_RPM", "missing_reads", "length_mismatches",
]

SCANNER = None
WINDOW = 80


def init_worker(barcodes: Sequence[dict], max_sub: int, max_ins: int, max_del: int, window: int) -> None:
    global SCANNER, WINDOW
    SCANNER = OccurrenceScanner(list(barcodes), max_sub, max_ins, max_del)
    WINDOW = window


def rescan_task(task: Tuple[str, Dict[str, List[Tuple[str, int, bool]]]]) -> dict:
    """Stream one FASTQ and rescan only the selected reads.

    selected maps read_id -> [(tool, db_read_length, legacy_internal), ...].
    """
    path, selected = task
    counts: Dict[str, Counter] = defaultdict(Counter)
    found = set()
    reclassified = []
    for header, seq, _, _ in fastq_reader(Path(path)):
        read_id = extract_read_id(header)
        entries = selected.get(read_id)
        if entries is None or read_id in found:
            continue
        found.add(read_id)
        internal = internal_count(SCANNER.scan(seq), len(seq), WINDOW) > 0
        for tool, db_length, legacy_internal in entries:
            tally = counts[tool]
            if db_length != len(seq):
                tally["length_mismatches"] += 1
            if legacy_internal:
                tally["internal_reads_found"] += 1
                if not internal:
                    tally["internal_reads_lost"] += 1
                    reclassified.append((tool, read_id, "internal_lost"))
            else:
                tally["reads_found"] += 1
                if internal:
                    tally["newly_internal_reads"] += 1
                    reclassified.append((tool, read_id, "newly_internal"))
    missing = [(tool, read_id) for read_id, entries in selected.items() if read_id not in found
               for tool, _, _ in entries]
    return {"path": path, "counts": {tool: dict(tally) for tool, tally in counts.items()},
            "missing": missing, "reclassified": reclassified}


def remap(path: str, path_map: Sequence[Tuple[str, str]]) -> str:
    for old, new in path_map:
        if path.startswith(old):
            return new + path[len(old):]
    return path


def open_readonly(db: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{db.resolve()}?mode=ro", uri=True)


def load_selection(db: Path, path_map, tools=None):
    """Return totals, legacy internal counts and {source_file: {read_id: [(tool, len, internal)]}}."""
    con = open_readonly(db)
    totals: Counter = Counter()
    legacy_internal: Counter = Counter()
    for tool, pattern_class, n in con.execute(
            "SELECT tool, pattern_class, COUNT(*) FROM reads GROUP BY tool, pattern_class"):
        if tools and tool not in tools:
            continue
        totals[tool] += n
        if pattern_class.endswith("__internal"):
            legacy_internal[tool] += n
    selection: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
    rescanned: Counter = Counter()
    rechecked: Counter = Counter()
    query = ("SELECT tool, source_file, read_id, read_length, pattern_class FROM reads "
             "WHERE pattern_class NOT IN (?, ?)")
    for tool, source_file, read_id, read_length, pattern_class in con.execute(query, NOT_RESCANNED):
        if tools and tool not in tools:
            continue
        is_internal = pattern_class.endswith("__internal")
        (rechecked if is_internal else rescanned)[tool] += 1
        selection[remap(source_file, path_map)][read_id].append((tool, read_length, is_internal))
    con.close()
    return totals, legacy_internal, rescanned, rechecked, selection


def split_tasks(selection, reads_per_task: int):
    """Split large files by a stable hash of the read ID, so one merged FASTQ still uses many cores.

    Each subtask re-reads the whole file (served from the page cache) but scans only its share.
    """
    tasks = []
    for path, reads in selection.items():
        parts = max(1, -(-len(reads) // reads_per_task))
        if parts == 1:
            tasks.append((path, dict(reads)))
            continue
        shares = [dict() for _ in range(parts)]
        for read_id, entries in reads.items():
            shares[zlib.crc32(read_id.encode()) % parts][read_id] = entries
        tasks.extend((path, share) for share in shares if share)
    return tasks


def rpm(n: int, total: int) -> float:
    return 1e6 * n / total if total else float("nan")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, required=True, help="comparison.sqlite of a best_hit run (opened read-only)")
    ap.add_argument("--barcode-csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True, help="Output CSV (one row per tool)")
    ap.add_argument("--experiment", default=None,
                    help="Experiment label for the CSV (default: database grandparent folder without _demux)")
    ap.add_argument("--terminal-window", type=int, default=80)
    ap.add_argument("--max-substitutions", type=int, default=1)
    ap.add_argument("--max-insertions", type=int, default=1)
    ap.add_argument("--max-deletions", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--reads-per-task", type=int, default=10000,
                    help="Selected reads per work item; larger files are split by read-ID hash")
    ap.add_argument("--tools", nargs="+", default=None, help="Only these tools/conditions (default: all)")
    ap.add_argument("--path-map", action="append", default=[], metavar="OLD=NEW",
                    help="Rewrite a source_file prefix, e.g. when the database was built on another host")
    ap.add_argument("--reclassified-out", type=Path, default=None,
                    help="Optional LOCAL list of reclassified read IDs (contains read identifiers)")
    args = ap.parse_args(argv)

    path_map = []
    for item in args.path_map:
        if "=" not in item:
            ap.error(f"--path-map needs OLD=NEW, got {item!r}")
        path_map.append(tuple(item.split("=", 1)))
    experiment = args.experiment or args.db.resolve().parent.parent.name.removesuffix("_demux")

    t0 = time.time()
    barcodes = load_barcodes(args.barcode_csv)
    totals, legacy_internal, rescanned, rechecked, selection = load_selection(
        args.db, path_map, set(args.tools) if args.tools else None)
    missing_files = [path for path in selection if not Path(path).exists()]
    if missing_files:
        raise SystemExit(f"{len(missing_files)} source file(s) not found, e.g. {missing_files[0]}; "
                         "use --path-map OLD=NEW")
    tasks = split_tasks(selection, args.reads_per_task)
    print(f"[{experiment}] selected {sum(rescanned.values()):,} non-internal and "
          f"{sum(rechecked.values()):,} internal reads in {len(selection)} files "
          f"({len(tasks)} work items, {time.time() - t0:.0f}s)", flush=True)

    counts: Dict[str, Counter] = defaultdict(Counter)
    missing: List[Tuple[str, str]] = []
    reclassified: List[Tuple[str, str, str]] = []
    with ProcessPoolExecutor(max_workers=args.threads, initializer=init_worker,
                             initargs=(barcodes, args.max_substitutions, args.max_insertions,
                                       args.max_deletions, args.terminal_window)) as ex:
        futures = [ex.submit(rescan_task, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            result = fut.result()
            for tool, tally in result["counts"].items():
                counts[tool].update(tally)
            missing.extend(result["missing"])
            reclassified.extend(result["reclassified"])
            if i % 10 == 0 or i == len(futures):
                print(f"[{experiment}] work items done: {i}/{len(futures)} ({time.time() - t0:.0f}s)", flush=True)

    missing_by_tool = Counter(tool for tool, _ in missing)
    for tool, n in sorted(missing_by_tool.items()):
        print(f"WARNING [{experiment}/{tool}] {n:,} selected read(s) not found in their source FASTQ",
              file=sys.stderr)
    for tool in sorted(counts):
        if counts[tool]["length_mismatches"]:
            print(f"WARNING [{experiment}/{tool}] {counts[tool]['length_mismatches']:,} read(s) "
                  "differ in length from the database", file=sys.stderr)

    rows = []
    for tool in sorted(totals):
        c = counts[tool]
        corrected = legacy_internal[tool] + c["newly_internal_reads"] - c["internal_reads_lost"]
        rows.append({
            "experiment": experiment,
            "tool": tool,
            "total_reads": totals[tool],
            "legacy_internal_reads": legacy_internal[tool],
            "legacy_internal_RPM": rpm(legacy_internal[tool], totals[tool]),
            "rescanned_reads": rescanned[tool],
            "reads_found": c["reads_found"],
            "newly_internal_reads": c["newly_internal_reads"],
            "fraction_of_rescanned_reclassified": (c["newly_internal_reads"] / rescanned[tool]
                                                   if rescanned[tool] else 0.0),
            "internal_reads_rechecked": rechecked[tool],
            "internal_reads_found": c["internal_reads_found"],
            "internal_reads_lost": c["internal_reads_lost"],
            "corrected_internal_reads": corrected,
            "corrected_internal_RPM": rpm(corrected, totals[tool]),
            "missing_reads": missing_by_tool[tool],
            "length_mismatches": c["length_mismatches"],
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    if args.reclassified_out:
        args.reclassified_out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.reclassified_out, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["tool", "read_id", "change"])
            writer.writerows(sorted(reclassified))

    print(f"\n{'tool':6s} {'legacy_RPM':>11s} {'newly_int':>10s} {'lost':>6s} {'corrected_RPM':>14s}")
    for row in rows:
        print(f"{row['tool']:6s} {row['legacy_internal_RPM']:11.2f} {row['newly_internal_reads']:10,d} "
              f"{row['internal_reads_lost']:6,d} {row['corrected_internal_RPM']:14.2f}")
    print(f"Written: {args.out} ({time.time() - t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
