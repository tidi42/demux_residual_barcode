#!/usr/bin/env python3
"""Sensitivity of the genomic-element composition to the fixed keyword precedence.

Count titles matching multiple category keywords and compare category composition
when plasmid/rRNA precedence is reversed.
"""
import argparse
import re
import pandas as pd

COLS = ["qseqid", "saccver", "staxids", "sscinames", "pident", "length", "qlen",
        "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "stitle"]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--barcode-hits", required=True)
parser.add_argument("--motif-hits", required=True)
arguments = parser.parse_args()

RULES = [
    ("Plasmid",        r"plasmid"),
    ("rRNA gene",      r"16s|18s|ribosomal rna|rrna"),
    ("Phage / virus",  r"phage|virus|viral"),
    ("Organelle",      r"mitochondri|chloroplast|plastid|apicoplast"),
    ("Bacterial chr.", r"chromosome|complete genome|genome assembly"),
]


def classify(title, order):
    s = str(title).lower()
    for name in order:
        if re.search(dict(RULES)[name], s):
            return name
    return "Other/uncultured"


def load(panel, fname):
    df = pd.read_csv(fname, sep="\t", names=COLS, header=None)
    df = df[(df.mismatch == 0) & (df.pident == 100) & (df.qcovs == 100)].copy()
    df["panel"] = panel
    return df


perfect = pd.concat([load("Barcode only", arguments.barcode_hits),
                     load("Adapter + barcode", arguments.motif_hits)],
                    ignore_index=True)
bc = perfect[perfect.panel == "Barcode only"]

print(f"perfect HSPs: total {len(perfect)}, barcode-only {len(bc)}")

# --- how many titles match more than one category keyword -------------------
hits = pd.DataFrame({name: bc.stitle.str.lower().str.contains(pat, regex=True)
                     for name, pat in RULES})
n_match = hits.sum(axis=1)
print("\nnumber of category keywords matched per subject title (barcode-only HSPs):")
for k, v in n_match.value_counts().sort_index().items():
    print(f"  {k} keyword(s): {v:6d} HSPs ({100*v/len(bc):5.1f}%)")

multi = hits[n_match >= 2]
print(f"\nambiguous HSPs (>= 2 keywords): {len(multi)} ({100*len(multi)/len(bc):.1f}%)")
print("co-occurring keyword pairs:")
combos = multi.apply(lambda r: " + ".join(sorted(multi.columns[r.values])), axis=1)
for combo, n in combos.value_counts().items():
    print(f"  {combo:45s} {n:5d}")

print(f"\nplasmid AND rRNA in the same title: "
      f"{int((hits['Plasmid'] & hits['rRNA gene']).sum())} HSPs")

# --- composition under the two precedence orders ----------------------------
default = [n for n, _ in RULES]
swapped = ["rRNA gene", "Plasmid"] + default[2:]

comp = {}
for label, order in [("default (plasmid first)", default), ("swapped (rRNA first)", swapped)]:
    lab = bc.stitle.map(lambda t: classify(t, order))
    comp[label] = 100 * lab.value_counts() / len(bc)

out = pd.DataFrame(comp).fillna(0.0)
out["delta (pp)"] = out.iloc[:, 1] - out.iloc[:, 0]
print("\ngenomic-element composition, barcode-only perfect HSPs (%):")
print(out.round(2).sort_values(out.columns[0], ascending=False).to_string())
