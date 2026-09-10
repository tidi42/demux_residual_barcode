"""Extract per-read residual-barcode multiplicity for a GraphPad histogram.

The manuscript reports "% of reads barcode-positive", which collapses a read
carrying one barcode and a read carrying five into the same bucket. This script
recovers the multiplicity distribution from the read-level SQLite databases.

Two views are produced:
  * residual_count           - how many barcode matches sit on a read
  * residual_unique_bc_count - how many DISTINCT barcode IDs sit on a read
                               (>1 is the crosstalk/chimera-relevant case)
"""
import argparse
import csv
import os
import sqlite3
import statistics as st
import time

parser = argparse.ArgumentParser(description="Summarise retained barcode-pattern multiplicity.")
parser.add_argument("--base-dir", default="outputs")
parser.add_argument("--experiments", nargs="+", required=True)
parser.add_argument("--out", default="figures/summary_figures/graphpad")
arguments = parser.parse_args()
EXPS = arguments.experiments
CONDS = ["c0", "c1", "c2", "c3", "c4", "c5", "c6"]
LABEL = {"c0": "C0", "c1": "C1", "c2": "C2", "c3": "C3",
         "c4": "C4", "c5": "C5", "c6": "C6"}
CAP_N, CAP_U = 5, 3          # top bins are "5 or more" / "3 or more"
OUT = arguments.out

# raw[exp][cond][n_matches] = reads ; uniq[exp][cond][n_distinct] = reads
raw, uniq, totals = {}, {}, {}

for e in EXPS:
    db = os.path.join(arguments.base_dir, f"{e}_demux", "cmp_all", "comparison.sqlite")
    if not os.path.exists(db):
        print(f"  !! missing {db}")
        continue
    t0 = time.time()
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    raw[e], uniq[e], totals[e] = {}, {}, {}
    q = (f"SELECT tool, MIN(residual_count,{CAP_N}), "
         f"MIN(residual_unique_bc_count,{CAP_U}), COUNT(*) "
         "FROM reads GROUP BY 1,2,3")
    for tool, n, u, cnt in con.execute(q):
        if tool not in CONDS:
            continue
        raw[e].setdefault(tool, {}).setdefault(n, 0)
        raw[e][tool][n] += cnt
        totals[e][tool] = totals[e].get(tool, 0) + cnt
        if n > 0:                       # unique-ID view only for positive reads
            uniq[e].setdefault(tool, {}).setdefault(u, 0)
            uniq[e][tool][u] += cnt
    con.close()
    print(f"  {e}: {time.time()-t0:.0f}s")

os.makedirs(OUT, exist_ok=True)
BINS_N = list(range(1, CAP_N + 1))       # 1..4, then "5+"
BINS_U = list(range(1, CAP_U + 1))       # 1, 2, then "3+"


def nlabel(n, cap):
    return f"{n}" if n < cap else f"{cap}+"


def pct_of_positive(d, bins, cap):
    """Percent of BARCODE-POSITIVE reads falling in each multiplicity bin."""
    pos = sum(v for k, v in d.items() if k > 0)
    if not pos:
        return {b: float("nan") for b in bins}
    return {b: d.get(b, 0) / pos * 100 for b in bins}


# ---- 1. per-experiment values, GraphPad "Grouped" layout ------------------
with open(f"{OUT}/graphpad_figH1_barcode_multiplicity_per_experiment.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Barcodes per read"] +
               [f"{LABEL[c]}_{e}" for c in CONDS for e in EXPS])
    for b in BINS_N:
        row = [nlabel(b, CAP_N)]
        for c in CONDS:
            for e in EXPS:
                row.append(round(pct_of_positive(raw[e].get(c, {}), BINS_N, CAP_N)[b], 4))
        w.writerow(row)

# ---- 2. mean composition, GraphPad stacked-bar layout ---------------------
with open(f"{OUT}/graphpad_figH2_barcode_multiplicity_stacked.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Condition"] + [f"{nlabel(b, CAP_N)} barcode(s)" for b in BINS_N])
    for c in CONDS:
        vals = [pct_of_positive(raw[e].get(c, {}), BINS_N, CAP_N) for e in EXPS]
        w.writerow([LABEL[c]] +
                   [round(st.mean(v[b] for v in vals), 3) for b in BINS_N])

# ---- 3. distinct barcode IDs per positive read ---------------------------
with open(f"{OUT}/graphpad_figH3_distinct_barcode_ids_per_read.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Distinct barcode IDs per positive read"] +
               [f"{LABEL[c]}_{e}" for c in CONDS for e in EXPS])
    for b in BINS_U:
        row = [nlabel(b, CAP_U)]
        for c in CONDS:
            for e in EXPS:
                d = uniq[e].get(c, {})
                tot = sum(d.values())
                row.append(round(d.get(b, 0) / tot * 100, 4) if tot else "")
        w.writerow(row)

# ---- 4. summary the reviewer can quote -----------------------------------
with open(f"{OUT}/summary_barcode_multiplicity.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["condition", "n_exp", "pct_reads_positive_mean",
                "pct_of_positive_with_1_mean", "pct_of_positive_with_1_sd",
                "pct_of_positive_with_2_mean", "pct_of_positive_with_2_sd",
                "pct_of_positive_with_3plus_mean",
                "mean_barcodes_per_positive_read", "median_barcodes_per_positive_read",
                "pct_of_positive_multi_distinct_ids_mean"])
    for c in CONDS:
        p1, p2, p3, posr, means, meds, multi = [], [], [], [], [], [], []
        for e in EXPS:
            d = raw[e].get(c, {})
            if not d:
                continue
            pos = sum(v for k, v in d.items() if k > 0)
            tot = totals[e].get(c, 0)
            if not pos or not tot:
                continue
            pc = pct_of_positive(d, BINS_N, CAP_N)
            p1.append(pc[1]); p2.append(pc[2])
            p3.append(sum(pc[b] for b in BINS_N if b >= 3))
            posr.append(pos / tot * 100)
            means.append(sum(k * v for k, v in d.items() if k > 0) / pos)
            # median of a counted distribution
            run, half, med = 0, pos / 2, None
            for k in sorted(b for b in d if b > 0):
                run += d[k]
                if run >= half:
                    med = k
                    break
            meds.append(med)
            u = uniq[e].get(c, {})
            ut = sum(u.values())
            multi.append(sum(v for k, v in u.items() if k >= 2) / ut * 100 if ut else 0)
        w.writerow([LABEL[c], len(p1), round(st.mean(posr), 4),
                    round(st.mean(p1), 3), round(st.stdev(p1), 3),
                    round(st.mean(p2), 3), round(st.stdev(p2), 3),
                    round(st.mean(p3), 3),
                    round(st.mean(means), 4), st.median(meds),
                    round(st.mean(multi), 3)])

print("\nwritten to", OUT)
for f in sorted(os.listdir(OUT)):
    if "multiplicity" in f or "distinct_barcode" in f:
        print("  ", f)
