#!/usr/bin/env python3
"""Release only cross-experiment and synthetic summaries of a detector pilot."""

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(source, output):
    output.mkdir(parents=True, exist_ok=True)
    rows = read_rows(source / "detector_sensitivity_summary_LOCAL.csv")
    groups = defaultdict(list)
    for row in rows:
        if row["group"] != "synthetic":
            groups[row["condition"], row["detector"], row["terminal_window"], row["length_stratum"]].append(row)
    totals = []
    for key, group in sorted(groups.items()):
        totals.append(dict(condition=key[0], detector=key[1], terminal_window=int(key[2]), length_stratum=key[3],
                           experiments_represented=len(group),
                           sampled_reads=sum(int(float(row["sample_reads"])) for row in group),
                           sampled_positive_reads=sum(int(float(row["positive_reads"])) for row in group),
                           sampled_internal_reads=sum(int(float(row["internal_reads"])) for row in group),
                           sampled_retained_matches=sum(int(float(row["retained_matches"])) for row in group),
                           mean_experiment_weighted_internal_RPM=statistics.mean(
                               float(row["weighted_internal_RPM"]) for row in group)))
    if totals:
        write_rows(output / "study_pilot_aggregate.csv", totals)
        paired = []
        windows = sorted({row["terminal_window"] for row in rows}, key=int)
        for first, second in [("c2", "c4"), ("c3", "c5")]:
            for detector in ("legacy_best_match", "occurrences_nonoverlap"):
                for window in windows:
                    for stratum in ("all_below_5000", "<250", "250-499", "500-999", "1000-4999"):
                        values = {}
                        for condition in (first, second):
                            values[condition] = {row["group"]: float(row["weighted_internal_RPM"]) for row in rows
                                                 if row["condition"] == condition and row["detector"] == detector
                                                 and row["terminal_window"] == window and row["length_stratum"] == stratum}
                        common = sorted(values[first].keys() & values[second].keys())
                        differences = [values[first][group] - values[second][group] for group in common]
                        paired.append(dict(comparison=first+"_vs_"+second, detector=detector,
                                           terminal_window=int(window), length_stratum=stratum,
                                           paired_experiments=len(common),
                                           mean_paired_difference_RPM=statistics.mean(differences) if differences else "",
                                           experiments_first_higher=sum(value > 0 for value in differences),
                                           experiments_equal=sum(value == 0 for value in differences),
                                           experiments_second_higher=sum(value < 0 for value in differences)))
        write_rows(output / "paired_pilot_aggregate.csv", paired)
    recovery = source / "synthetic_recovery.csv"
    if recovery.exists():
        write_rows(output / "synthetic_recovery.csv", read_rows(recovery))
    background = [{column: row[column] for column in ["condition", "detector", "terminal_window", "length_stratum",
                                                      "sample_reads", "positive_reads", "internal_reads", "retained_matches"]}
                  for row in rows if row["group"] == "synthetic" and row["condition"] == "barcode_free"]
    if background:
        write_rows(output / "synthetic_background.csv", background)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    arguments = parser.parse_args()
    summarize(arguments.source, arguments.out)
    print("Wrote aggregate summaries without experiment identifiers or replicate columns")


if __name__ == "__main__":
    main()