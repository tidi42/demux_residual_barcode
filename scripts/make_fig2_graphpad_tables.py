#!/usr/bin/env python3
"""GraphPad import tables for Figure 2, built from the COMPLETE BLAST output.

Figure 2C is drawn in Prism, so it needs a table the author can paste in
directly: rows are barcodes, columns are the six genomic-element categories.
Both a percentage table (100% stacked bars) and a count table (absolute
stacked bars) are written, ordered identically to Figure 2B so the two panels
line up.

Also writes a species-level aggregation of the Figure 2A data. The published
"organism" is the subject description up to the first comma, which is
strain-level: six Photobacterium leiognathi strains are six separate bubbles.
Collapsing to species changes what the panel shows, so both are provided and
the choice is left to the author.

Supply both BLAST tables explicitly; no input or result data are bundled.
"""
import argparse
import re
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser(description="Aggregate the operational barcode-database screen.")
parser.add_argument("--barcode-hits", type=Path, required=True)
parser.add_argument("--motif-hits", type=Path, required=True)
parser.add_argument("--out", type=Path, default=Path("figures/fig2"))
arguments = parser.parse_args()
RAW = {"Barcode only": arguments.barcode_hits, "Adapter + barcode": arguments.motif_hits}
OUT = arguments.out
OUT.mkdir(parents=True, exist_ok=True)
COLS = ["qseqid", "saccver", "staxids", "sscinames", "pident", "length", "qlen",
        "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "stitle"]
CATEGORIES = ["Bacterial chr.", "rRNA gene", "Plasmid", "Phage / virus",
              "Organelle", "Other/uncultured"]
RULES = [
    ("Plasmid", r"plasmid"),
    ("rRNA gene", r"16s|18s|ribosomal rna|rrna"),
    ("Phage / virus", r"phage|virus|viral"),
    ("Organelle", r"mitochondri|chloroplast|plastid|apicoplast"),
    ("Bacterial chr.", r"chromosome|complete genome|genome assembly"),
]


def classify_title(title):
    for category, pattern in RULES:
        if re.search(pattern, str(title), flags=re.IGNORECASE):
            return category
    return "Other/uncultured"


def load_perfect(path):
    df = pd.read_csv(path, sep="\t", names=COLS, header=None)
    df = df[(df.mismatch == 0) & (df.pident == 100) & (df.qcovs == 100)].copy()
    df["barcode"] = df.qseqid.str.extract(r"(NB\d+)")
    df["organism"] = df.stitle.str.replace(r",.*$", "", regex=True).str.strip()
    return df


QUALIFIER = {"clone", "strain", "isolate", "genome", "chromosome", "plasmid",
             "genotype", "serovar", "voucher"}


def species_of(name):
    """Genus + species, stopping before any isolate identifier."""
    unc = bool(re.match(r"^[Uu]ncultured\s+", name))
    tok = re.sub(r"^[Uu]ncultured\s+", "", name).split()
    cut = next((i for i, t in enumerate(tok) if t.lower() in QUALIFIER), len(tok))
    tok = tok[:cut] or tok[:1]
    if len(tok) >= 3 and tok[1] in ("phage", "virus", "sp."):
        lab = " ".join(tok[:3])
    else:
        lab = " ".join(tok[:2])
    # "Uncultured bacterium clone X" would otherwise collapse to "bacterium"
    if unc and lab.lower() in ("bacterium", "archaeon", "organism", "prokaryote"):
        lab = "uncultured " + lab
    return lab


# ---------------------------------------------------------------- Figure 2C
element_frames = []
for panel, path in RAW.items():
    hits = load_perfect(path)
    hits["entry_type"] = hits.stitle.map(classify_title)
    elements = hits.groupby(["barcode", "entry_type"]).size().reset_index(name="n_hits")
    elements.insert(0, "panel", panel)
    elements["total_hits"] = elements.groupby("barcode").n_hits.transform("sum")
    elements["proportion"] = elements.n_hits / elements.total_hits
    element_frames.append(elements)
all_elements = pd.concat(element_frames, ignore_index=True)
all_elements.to_csv(f"{OUT}/Fig2C_genomic_element_FULLDATA.csv", index=False)
c = all_elements[all_elements.panel == "Barcode only"]
assert not c.empty, "No barcode-only HSPs passed the operational filter"
order = (c.groupby("barcode").total_hits.first()
         .sort_values(ascending=False).index.tolist())

counts = (c.pivot_table(index="barcode", columns="entry_type", values="n_hits",
                        aggfunc="sum", fill_value=0)
          .reindex(index=order, columns=CATEGORIES, fill_value=0)
          .astype(int))
pct = (100 * counts.T / counts.sum(axis=1)).T.round(4)

assert counts.values.sum() == c.n_hits.sum(), counts.values.sum()
assert pct.sum(axis=1).round(3).eq(100).all()

counts.rename_axis("Barcode").to_csv(f"{OUT}/graphpad_Fig2C_element_counts.csv")
pct.rename_axis("Barcode").to_csv(f"{OUT}/graphpad_Fig2C_element_percent.csv")

# ---------------------------------------------------------------- Figure 2A
frames = []
subject_frames = []
for panel, path in RAW.items():
    q = load_perfect(path)
    subjects = (q.groupby("organism")
                .agg(breadth_distinct_barcodes=("barcode", "nunique"),
                     depth_hsp=("saccver", "size"),
                     depth_distinct_subjects=("saccver", "nunique"),
                     max_bitscore=("bitscore", "max"))
                .reset_index())
    subjects.insert(0, "panel", panel)
    subjects["risk_score_hsp"] = subjects.breadth_distinct_barcodes * subjects.depth_hsp
    subject_frames.append(subjects)
    q["species"] = q.organism.map(species_of)
    agg = (q.groupby("species")
           .agg(breadth_distinct_barcodes=("barcode", "nunique"),
                depth_hsp=("saccver", "size"),
                depth_distinct_subjects=("saccver", "nunique"),
                max_bitscore=("bitscore", "max"),
                n_strain_entries=("organism", "nunique"))
           .reset_index())
    agg.insert(0, "panel", panel)
    frames.append(agg)
    if panel == "Barcode only":
        p = q

pd.concat(subject_frames, ignore_index=True).to_csv(
    f"{OUT}/Fig2A_shared_vulnerability_FULLDATA.csv", index=False)
sp = pd.concat(frames, ignore_index=True)
sp["risk_score_hsp"] = sp.breadth_distinct_barcodes * sp.depth_hsp
sp.sort_values(["panel", "risk_score_hsp"], ascending=[True, False]).to_csv(
    f"{OUT}/Fig2A_shared_vulnerability_SPECIES.csv", index=False)
bo = sp[sp.panel == "Barcode only"]

# ------------------------------------------------- Figure 2B per-barcode table
# "N spp." must count species, not subject records: NB03 matches 31 records but
# only 17 distinct species.
ranks = []
for panel, path in RAW.items():
    q = load_perfect(path)
    q["species"] = q.organism.map(species_of)
    r = (q.groupby("barcode")
         .agg(n_perfect_hsp=("saccver", "size"),
              n_perfect_distinct_subjects=("saccver", "nunique"),
              n_organisms=("organism", "nunique"),
              n_species=("species", "nunique"))
         .reset_index())
    r.insert(0, "panel", panel)
    ranks.append(r)
rank = pd.concat(ranks, ignore_index=True)
rank["risk_class_by_distinct_subjects"] = pd.cut(
    rank.n_perfect_distinct_subjects, [-1, 9, 19, 10 ** 9],
    labels=["Low (1-9 hits)", "Moderate (10-19 hits)", "High (>=20 hits)"])
rank.sort_values(["panel", "n_perfect_distinct_subjects"],
                 ascending=[True, False]).to_csv(
    f"{OUT}/Fig2B_barcode_ranking_FULLDATA.csv", index=False)

print(f"Figure 2C: {len(counts)} barcodes x {len(CATEGORIES)} categories, "
      f"{counts.values.sum()} hits")
print(f"  {OUT}/graphpad_Fig2C_element_counts.csv")
print(f"  {OUT}/graphpad_Fig2C_element_percent.csv")
print("\noverall composition (% of retained HSPs):")
tot = 100 * counts.sum() / counts.values.sum()
for k, v in tot.items():
    print(f"  {k:18s} {v:5.2f}%")

print(f"\nFigure 2A species-level: {len(bo)} species "
      f"(vs {p.organism.nunique()} strain-level entries)")
print(f"  matched by >1 barcode: {(bo.breadth_distinct_barcodes >= 2).sum()} "
      f"(strain-level: {(p.groupby('organism').barcode.nunique() >= 2).sum()})")
print(f"  {OUT}/Fig2A_shared_vulnerability_SPECIES.csv")
print(bo.nlargest(8, "risk_score_hsp")[
    ["species", "breadth_distinct_barcodes", "depth_hsp", "risk_score_hsp"]
].to_string(index=False))

rb = rank[rank.panel == "Barcode only"]
print(f"\nFigure 2B table: {len(rb)} barcodes, {rb.n_perfect_hsp.sum()} HSPs")
print(f"  max species per barcode {rb.n_species.max()} "
      f"({rb.loc[rb.n_species.idxmax(), 'barcode']}); "
      f"max subject records {rb.n_organisms.max()}")
print(f"  {OUT}/Fig2B_barcode_ranking_FULLDATA.csv")
