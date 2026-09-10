#!/usr/bin/env python3
"""
Generate summary_f_log_terminal_vs_internal.csv — the log10-axis-ready
reformatting of Figure panel F (residual-barcode location rate).

Panel F values are computed EXACTLY as in Table S4:

    reads_per_million = class_reads / total_reads * 1e6

per residual location (terminal / internal / both_ends), per condition
(C0-C6), per experiment (n = 8). This script reads the verified Table S4
per-experiment source (figure_f_violin/supplement_f_per_experiment.csv) and,
for every (location, condition), emits:

  * the 8 per-experiment RPM values (exp_*), and
  * log10-appropriate summary statistics for a log10 y-axis:
        geometric mean (+ its log10), median (+ its log10), Q1, Q3, IQR,
        min, max.

Arithmetic mean +/- SD (used by the old linear stacked-bar panel F, and the
basis of the reviewer criticism) is deliberately omitted: it is not meaningful
on a log10 axis for a quantity that spans ~5 orders of magnitude across
conditions. Median/IQR and the geometric mean are the correct central-tendency
and spread measures on a log scale.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description="Calculate log-scale location-rate summaries.")
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--experiments", nargs="+", required=True)
parser.add_argument("--out", type=Path, required=True)
arguments = parser.parse_args()
SRC = arguments.input
OUT = arguments.out

GROUP = {
    "c0": "raw", "c1": "raw",
    "c2": "Dorado demux", "c3": "Dorado demux",
    "c4": "Barbell", "c5": "Barbell",
    "c6": "Dorado+Barbell",
}
CONDITIONS = ["c0", "c1", "c2", "c3", "c4", "c5", "c6"]
LOCATIONS = ["terminal", "internal", "both_ends"]
EXPERIMENTS = arguments.experiments


def r(x: float, nd: int = 4) -> float:
    return round(float(x), nd)


# (location, condition) -> {experiment: reads_per_million}
data: dict[tuple[str, str], dict[str, float]] = {}
with SRC.open() as fh:
    for row in csv.DictReader(fh):
        key = (row["location"], row["condition"])
        data.setdefault(key, {})[row["experiment"]] = float(row["reads_per_million"])

header = (
    ["location", "condition", "group", "stat_n",
     "geomean", "log10_geomean", "median", "log10_median",
     "q1", "q3", "iqr", "min", "max"]
    + [f"exp_{e}" for e in EXPERIMENTS]
)

out_rows: list[dict] = []
for loc in LOCATIONS:
    for cond in CONDITIONS:
        per_exp = data[(loc, cond)]
        vals = np.array([per_exp[e] for e in EXPERIMENTS], dtype=float)
        logs = np.log10(vals)  # all values > 0, so log10 is always defined
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        rec: dict = {
            "location": loc,
            "condition": cond,
            "group": GROUP[cond],
            "stat_n": int(vals.size),
            "geomean": r(10.0 ** logs.mean()),
            "log10_geomean": r(logs.mean()),
            "median": r(med),
            "log10_median": r(math.log10(med)),
            "q1": r(q1),
            "q3": r(q3),
            "iqr": r(q3 - q1),
            "min": r(vals.min()),
            "max": r(vals.max()),
        }
        for e in EXPERIMENTS:
            rec[f"exp_{e}"] = r(per_exp[e])
        out_rows.append(rec)

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=header)
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Wrote {OUT} ({len(out_rows)} rows: "
      f"{len(LOCATIONS)} locations x {len(CONDITIONS)} conditions)")
