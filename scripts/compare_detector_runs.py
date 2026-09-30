#!/usr/bin/env python3
"""Quick old-vs-new comparison of two residual-barcode result trees.

Compares a best_hit tree (default outputs/) with an occurrences rerun (default
outputs_occurrences/), both laid out as {base}/{experiment}_demux/cmp_all/, and writes a
Markdown report with:
  * provenance and a read-total sanity check (the same FASTQs must give the same totals)
  * per-condition location rates (internal / terminal / both-ends RPM, % barcode-positive)
  * pattern classes that only the occurrences detector can populate
  * paired Wilcoxon statistics before and after (supplement_f_paired_stats.csv)
  * the targeted rescan cross-check (rescan_masked_internal.py vs the full rerun)
  * the decoy background, where decoy_background.csv exists
  * Fig. 1D composition (pattern view) and the occurrence view
  * which GraphPad CSVs changed, with the largest absolute cell change

The report contains aggregate numbers only (no read identifiers), but it describes
study results: keep it with the local outputs, not in Git.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

EXPERIMENTS = []
CONDITIONS = ["c0", "c1", "c2", "c3", "c4", "c5", "c6"]
LOCATIONS = ["internal", "both_ends", "terminal"]


def read_csv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def location(pattern_class: str) -> str:
    for suffix in LOCATIONS:
        if pattern_class.endswith("__" + suffix):
            return suffix
    return "none"


def load_tree(base: Path, subdir: str, experiments):
    """data[exp][cond] = {total, positive_pct, internal, both_ends, terminal, classes{}}."""
    data = {}
    for e in experiments:
        cmp_dir = base / f"{e}_demux" / subdir
        if not (cmp_dir / "tool_summary.csv").exists():
            continue
        data[e] = {}
        for r in read_csv(cmp_dir / "tool_summary.csv"):
            data[e][r["tool"]] = {"total": int(r["total_reads"]),
                                  "positive_pct": float(r["percent_with_residual_barcode"]),
                                  "internal": 0, "both_ends": 0, "terminal": 0, "classes": {}}
        for r in read_csv(cmp_dir / "residual_pattern_summary.csv"):
            d = data[e].get(r["tool"])
            if d is None:
                continue
            n = int(r["reads"])
            d["classes"][r["pattern_class"]] = n
            loc = location(r["pattern_class"])
            if loc != "none":
                d[loc] += n
    return data


def rpm(d, key):
    return 1e6 * d[key] / d["total"] if d["total"] else float("nan")


def geo_mean(values):
    values = [v for v in values if v and v > 0 and not math.isnan(v)]
    return math.exp(st.mean(math.log(v) for v in values)) if values else float("nan")


def fmt(value, digits=2):
    if value == "" or value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return "n/a"
    if isinstance(value, int):
        return f"{value:,}"
    return f"{value:,.{digits}f}"


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def compare_csv(before: Path, after: Path):
    """Positional cell comparison of two CSVs -> (changed, compared, max_change, notes).

    Cells are matched by (row, column) position, so a renamed header does not count as a
    change; columns or rows present in only one file are reported in notes.
    """
    a = list(csv.reader(open(before, newline="")))
    b = list(csv.reader(open(after, newline="")))
    notes = []
    renamed = [i for i, (x, y) in enumerate(zip(a[0], b[0])) if x != y]
    if renamed:
        notes.append(f"{len(renamed)} header label(s) renamed")
    if len(b[0]) != len(a[0]):
        notes.append(f"{len(b[0]) - len(a[0]):+d} column(s)")
    if len(b) != len(a):
        notes.append(f"{len(b) - len(a):+d} row(s)")
    changed = compared = 0
    biggest = 0.0
    for row_a, row_b in zip(a[1:], b[1:]):
        for x, y in zip(row_a, row_b):
            if _is_number(x) and _is_number(y):
                compared += 1
                delta = abs(float(x) - float(y))
                if delta > 1e-9:
                    changed += 1
                    biggest = max(biggest, delta)
            elif x != y:
                compared += 1
                changed += 1
    return changed, compared, biggest, notes


def _is_number(cell):
    try:
        float(cell)
        return True
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-base", type=Path, default=Path("outputs"))
    ap.add_argument("--new-base", type=Path, default=Path("outputs_occurrences"))
    ap.add_argument("--cmp-subdir", default="cmp_all")
    ap.add_argument("--experiments", nargs="+", required=True)
    ap.add_argument("--rescan-dir", type=Path, default=Path("outputs_occurrences/rescan"))
    ap.add_argument("--old-figures", type=Path, default=Path("figures/summary_figures"))
    ap.add_argument("--new-figures", type=Path, default=Path("figures/summary_figures_occurrences"))
    ap.add_argument("--out", type=Path, default=Path("figures/summary_figures_occurrences/detector_comparison_report.md"))
    args = ap.parse_args()

    old = load_tree(args.old_base, args.cmp_subdir, args.experiments)
    new = load_tree(args.new_base, args.cmp_subdir, args.experiments)
    exps = [e for e in args.experiments if e in old and e in new]
    missing = [e for e in args.experiments if e not in exps]
    out = []
    out.append("# Residual-barcode detector comparison: best_hit (old) vs occurrences (new)\n\n")
    out.append(f"Old tree: `{args.old_base}/{{exp}}_demux/{args.cmp_subdir}` (best_hit). "
               f"New tree: `{args.new_base}/{{exp}}_demux/{args.cmp_subdir}`.\n\n")
    out.append(f"Experiments compared: {', '.join(exps) or 'none'}")
    out.append(f"; not yet available: {', '.join(missing)}.\n\n" if missing else ".\n\n")

    # ── Provenance ────────────────────────────────────────────────────────────
    provenance = [json.loads((args.new_base / f"{e}_demux" / args.cmp_subdir / "detector_provenance.json").read_text())
                  for e in exps if (args.new_base / f"{e}_demux" / args.cmp_subdir / "detector_provenance.json").exists()]
    if provenance:
        p = provenance[0]
        same = all(q["scanner_sha256"] == p["scanner_sha256"] and q["production_script_sha256"] ==
                   p["production_script_sha256"] for q in provenance)
        out.append("## Provenance of the new results\n\n")
        out.append(f"- Detector: `{p['detector']}`; edit budget {p['edit_budget']}; terminal window "
                   f"{p['terminal_window']} bp; regex {p['regex_version']}\n")
        out.append(f"- Scanner `{p['scanner_source']}` SHA-256 `{p['scanner_sha256'][:16]}…`; production script "
                   f"SHA-256 `{p['production_script_sha256'][:16]}…` "
                   f"({'identical for all experiments' if same else 'DIFFERS between experiments'})\n")
        out.append(f"- `residual_count` now counts {p['residual_count_meaning']}\n\n")

    # ── Sanity: totals ────────────────────────────────────────────────────────
    mismatches = [(e, c, old[e][c]["total"], new[e].get(c, {}).get("total"))
                  for e in exps for c in CONDITIONS if c in old[e] and old[e][c]["total"] != new[e].get(c, {}).get("total")]
    out.append("## Sanity check\n\n")
    if mismatches:
        out.append("**Read totals differ between the trees** (the inputs are not identical):\n\n")
        out.append(table(["experiment", "condition", "old total", "new total"],
                         [(e, c.upper(), fmt(a), fmt(b)) for e, c, a, b in mismatches]))
    else:
        n = sum(old[e][c]["total"] for e in exps for c in CONDITIONS if c in old[e])
        out.append(f"Total reads are identical in every experiment and condition ({n:,} reads).\n")
    positive_changed = [(e, c) for e in exps for c in CONDITIONS if c in old[e] and
                        abs(old[e][c]["positive_pct"] - new[e][c]["positive_pct"]) > 1e-9]
    out.append("Barcode-positive read fractions (Fig. 1E) are "
               + ("identical in every experiment and condition, as expected: both detectors agree on "
                  "whether a read has any match.\n\n" if not positive_changed else
                  f"**different in {len(positive_changed)} experiment/condition pairs**, which is unexpected.\n\n"))

    # ── Location rates ────────────────────────────────────────────────────────
    out.append("## Location rates per condition (geometric mean over experiments, RPM)\n\n")
    rows = []
    for c in CONDITIONS:
        row = [c.upper()]
        for loc in LOCATIONS:
            a = geo_mean([rpm(old[e][c], loc) for e in exps if c in old[e]])
            b = geo_mean([rpm(new[e][c], loc) for e in exps if c in new[e]])
            row += [fmt(a), fmt(b), fmt(b / a if a else float("nan"), 3)]
        rows.append(row)
    out.append(table(["condition"] + [f"{loc} {x}" for loc in LOCATIONS for x in ("old", "new", "new/old")], rows))
    out.append("\nReads move mainly from terminal/both-ends to internal (a hidden internal copy is found) "
               "or from terminal to both-ends (a hidden copy at the other end). Very rarely the greedy overlap "
               "rule moves an internal best-hit match across the terminal-window boundary (see "
               "`internal_reads_lost` below). Total barcode-positive reads do not change.\n\n")

    out.append("### Internal RPM per experiment (primary endpoint)\n\n")
    rows = []
    for e in exps:
        row = [e]
        for c in ("c2", "c3", "c4", "c5"):
            a, b = rpm(old[e][c], "internal"), rpm(new[e][c], "internal")
            row += [f"{fmt(a)} → {fmt(b)}"]
        rows.append(row)
    out.append(table(["experiment", "C2", "C3", "C4", "C5"], rows))
    out.append("\n")

    # ── New classes ───────────────────────────────────────────────────────────
    out.append("## Pattern classes only the occurrences detector can populate\n\n")
    rows = []
    for c in CONDITIONS:
        counts = defaultdict(int)
        for e in exps:
            for pc, n in new[e].get(c, {}).get("classes", {}).items():
                if pc.startswith("multiple_same_barcode_same_direction"):
                    counts[pc.split("__")[1]] += n
        total = sum(new[e][c]["total"] for e in exps if c in new[e])
        rows.append([c.upper()] + [fmt(counts[loc]) for loc in LOCATIONS] +
                    [fmt(1e6 * sum(counts.values()) / total if total else float("nan"))])
    out.append(table(["condition", "…same_direction__internal", "…__both_ends", "…__terminal", "all (RPM, pooled)"], rows))
    out.append("\nReads summed over experiments. These reads carry repeated copies of one barcode "
               "orientation, which best_hit collapsed to a single match.\n\n")

    # ── Paired statistics ─────────────────────────────────────────────────────
    old_paired = args.old_figures / "figure_f_violin" / "supplement_f_paired_stats.csv"
    new_paired = args.new_figures / "figure_f_violin" / "supplement_f_paired_stats.csv"
    if old_paired.exists() and new_paired.exists():
        out.append("## Paired Wilcoxon signed-rank tests (internal RPM, experiments paired)\n\n")
        a = {r["comparison"]: r for r in read_csv(old_paired)}
        b = {r["comparison"]: r for r in read_csv(new_paired)}
        rows = []
        for comp in a:
            if comp not in b:
                continue
            ra, rb = a[comp], b[comp]
            rows.append([comp, f"{ra['median_a_rpm']} / {ra['median_b_rpm']}", f"{rb['median_a_rpm']} / {rb['median_b_rpm']}",
                         f"{ra['median_fold_change']} ({ra['fold_change_min']}–{ra['fold_change_max']})",
                         f"{rb['median_fold_change']} ({rb['fold_change_min']}–{rb['fold_change_max']})",
                         f"{ra['n_experiments_a_gt_b']}/{ra['n']} → {rb['n_experiments_a_gt_b']}/{rb['n']}",
                         f"{ra['p_two_sided']} → {rb['p_two_sided']}"])
        out.append(table(["comparison", "median RPM A/B old", "median RPM A/B new", "fold old", "fold new",
                          "runs A>B", "p (two-sided)"], rows))
        out.append("\nWith eight paired experiments the smallest attainable two-sided p is 0.0078 "
                   "(all runs move in the same direction).\n\n")

    # ── Rescan cross-check ────────────────────────────────────────────────────
    rescans = sorted(args.rescan_dir.glob("*_rescan_masked_internal.csv")) if args.rescan_dir.exists() else []
    if rescans:
        out.append("## Targeted rescan cross-check\n\n")
        out.append("`rescan_masked_internal.py` corrected the old databases by rescanning only reads with a "
                   "non-internal best-hit match (and rechecking internal ones). Its corrected internal count must "
                   "equal the full occurrences rerun.\n\n")
        rows, agree, compared = [], 0, 0
        for path in rescans:
            for r in read_csv(path):
                e, c = r["experiment"], r["tool"]
                full = new.get(e, {}).get(c, {}).get("internal")
                match = "" if full is None else ("yes" if int(r["corrected_internal_reads"]) == full else "**no**")
                if full is not None:
                    compared += 1
                    agree += match == "yes"
                rows.append([e, c.upper(), fmt(int(r["legacy_internal_reads"])), fmt(int(r["newly_internal_reads"])),
                             fmt(int(r["internal_reads_lost"])), fmt(int(r["corrected_internal_reads"])),
                             fmt(full) if full is not None else "pending", match,
                             fmt(int(r["missing_reads"])), fmt(100 * float(r["fraction_of_rescanned_reclassified"]), 3)])
        out.append(f"Agreement with the full rerun: {agree}/{compared} experiment/condition pairs compared.\n\n")
        out.append(table(["exp", "cond", "legacy internal", "newly internal", "lost", "corrected", "full rerun",
                          "equal", "missing reads", "% of rescanned reclassified"], rows))
        out.append("\n")

    # ── Decoys ────────────────────────────────────────────────────────────────
    decoy_rows = []
    for e in exps:
        path = args.new_base / f"{e}_demux" / args.cmp_subdir / "decoy_background.csv"
        if path.exists():
            decoy_rows += [(e, r) for r in read_csv(path)]
    if decoy_rows:
        out.append("## Decoy background (shuffled barcode sets)\n\n")
        out.append("Each decoy set shuffles every barcode's bases (composition preserved; backward decoy = "
                   "reverse complement of the shuffled forward decoy) and scans the same reads with the same "
                   "detector. It measures the chance-match rate of random barcode-like sequences.\n\n")
        rows = [[e, r["tool"].upper(), r["decoy_set"], fmt(int(r["total_reads"])),
                 fmt(float(r["decoy_internal_RPM"])), fmt(float(r["real_internal_RPM"])),
                 fmt(float(r["decoy_terminal_RPM"])), fmt(float(r["decoy_both_ends_RPM"]))] for e, r in decoy_rows]
        out.append(table(["exp", "cond", "decoy set", "reads", "decoy internal RPM", "real internal RPM",
                          "decoy terminal RPM", "decoy both-ends RPM"], rows))
        out.append("\n")

    # ── Fig. 1D composition ───────────────────────────────────────────────────
    old_h2 = args.old_figures / "graphpad" / "graphpad_figH2_barcode_multiplicity_stacked.csv"
    new_h2 = args.new_figures / "graphpad" / "graphpad_figH2_barcode_multiplicity_stacked.csv"
    new_h2o = args.new_figures / "graphpad" / "graphpad_figH2_barcode_multiplicity_stacked_occurrences.csv"
    if old_h2.exists() and new_h2.exists():
        out.append("## Fig. 1D composition (% of positive reads, mean over experiments)\n\n")
        a, b = read_csv(old_h2), read_csv(new_h2)
        c = read_csv(new_h2o) if new_h2o.exists() else None
        ka, kb = list(a[0])[1:], list(b[0])[1:]
        rows = []
        for i, row in enumerate(a):
            cells = [row["Condition"]]
            for j in range(len(ka)):
                cell = f"{row[ka[j]]} → {b[i][kb[j]]}"
                if c:
                    cell += f" ({c[i][list(c[i])[1 + j]]})"
                cells.append(cell)
            rows.append(cells)
        out.append(table(["condition"] + [k.split()[0] for k in ka], rows))
        out.append("\nCells: old pattern view → new pattern view (new occurrence view). The pattern view counts "
                   "distinct barcode orientations and is the Fig. 1D quantity; the occurrence view counts copies.\n\n")

    # ── GraphPad files ────────────────────────────────────────────────────────
    old_gp, new_gp = args.old_figures / "graphpad", args.new_figures / "graphpad"
    if new_gp.exists():
        out.append("## GraphPad CSV files\n\n")
        rows = []
        for path in sorted(new_gp.glob("*.csv")):
            before = old_gp / path.name
            if not before.exists():
                rows.append([f"`{path.name}`", "new file", "", ""])
                continue
            changed, compared, biggest, notes = compare_csv(before, path)
            status = "unchanged" if not changed else f"{changed} of {compared} shared cells changed"
            rows.append([f"`{path.name}`", status, fmt(biggest, 4) if changed else "", "; ".join(notes)])
        out.append(table(["file", "status", "largest absolute change", "layout"], rows))
        out.append("\nCells are compared by position. \"Unchanged\" files need no new GraphPad import; "
                   "files marked new are additional views (occurrence counts).\n\n")

    out.append("## Unchanged by the detector\n\n")
    out.append("Figure 2 (barcode-reference BLAST audit, rendered in R), read yield, length distributions, N50 "
               "and Fig. 1E do not depend on the residual-barcode detector.\n")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(out))
    print(f"Report: {args.out}")


if __name__ == "__main__":
    main()
