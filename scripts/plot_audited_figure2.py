#!/usr/bin/env python3
"""Render Figure 2A/B from audited counts; retain composition as a table."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def load_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def style_axis(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=11)


def draw_taxa(axis, rows):
    selected = [row for row in rows if row["panel"] == "barcode_only"]
    axis.scatter([int(row["distinct_barcodes"]) for row in selected],
                 [int(row["distinct_accessions"]) for row in selected],
                 s=[8 * int(row["barcode_accession_pairs"]) for row in selected],
                 color="#0072B2", edgecolor="white", linewidth=0.6, alpha=0.65)
    for row in selected[:4]:
        horizontal = int(row["distinct_barcodes"])
        vertical = int(row["distinct_accessions"])
        axis.annotate(row["taxon_group"].replace("taxid:", "TaxID "), (horizontal, vertical),
                      xytext=(6, 7) if horizontal < 12 else (-4, -19),
                      textcoords="offset points", ha="left" if horizontal < 12 else "right",
                      fontsize=10)
    axis.set(xlim=(0.4, 15.5), ylim=(0.75, 105), yscale="log",
             xlabel="Distinct barcode IDs per TaxID group",
             ylabel="Distinct accessions per TaxID group")
    axis.set_xticks([1, 4, 7, 10, 14])
    axis.set_yticks([1, 2, 5, 10, 20, 50, 100])
    axis.yaxis.set_major_formatter(ScalarFormatter())
    axis.minorticks_off()
    axis.set_title("A  Exported TaxID groups", loc="left", fontweight="bold", fontsize=15, pad=16)
    axis.text(0, -0.28, "143 groups; coincident points overlap.\n"
              "Bubble area: barcode-accession pairs.\n"
              "Motif + barcode: no complete matches.", transform=axis.transAxes,
              fontsize=11, va="top", linespacing=1.5)
    style_axis(axis)


def draw_ranking(axis, rows):
    selected = [row for row in rows if row["panel"] == "barcode_only" and int(row["distinct_accessions"]) > 0]
    values = [int(row["distinct_accessions"]) for row in selected]
    positions = list(range(len(selected)))
    axis.barh(positions, values, color="#0072B2", height=0.76)
    axis.set_yticks(positions, [row["barcode"] for row in selected], fontsize=10)
    for position, value in zip(positions, values):
        axis.text(value + 0.9, position, str(value), va="center", fontsize=10)
    axis.invert_yaxis()
    axis.set_xlim(0, 81)
    axis.set_xticks([0, 20, 40, 60, 80])
    axis.set_xlabel("Distinct subject accessions per barcode", fontsize=11)
    axis.set_title("B  Barcode-accession counts", loc="left", fontweight="bold", fontsize=15, pad=16)
    axis.text(0, -0.062, "55/96 barcodes have complete matches; 41 have zero.\n"
              "Motif + barcode: zero for all 96 barcodes.", transform=axis.transAxes,
              fontsize=10, va="top", linespacing=1.5)
    axis.grid(axis="x", color="#dddddd", linewidth=0.5)
    axis.set_axisbelow(True)
    style_axis(axis)
    axis.tick_params(axis="y", length=0, labelsize=10)


def draw_composition(axis, rows):
    selected = [row for row in rows if row["panel"] == "barcode_only"]
    pair_rows = [row for row in selected if row["unit"] == "barcode_accession_pairs"]
    hsps = {row["category"]: row for row in selected if row["unit"] == "HSPs"}
    cells = [[row["category"].replace("Chromosome/genome", "Chromosome/\ngenome"),
              f"{row['count']}\n({float(row['percent']):.1f}%)",
              f"{hsps[row['category']]['count']}\n({float(hsps[row['category']]['percent']):.1f}%)"]
             for row in pair_rows]
    axis.axis("off")
    axis.set_title("C  Record-annotation composition", loc="left", fontweight="bold", fontsize=15, pad=12)
    table = axis.table(cellText=cells, colLabels=["Title category", "Pairs\n(n = 363)", "HSPs\n(n = 1,128)"],
                       cellLoc="left", colLoc="left", colWidths=[0.48, 0.26, 0.26], bbox=[0, 0.19, 1, 0.77])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_edgecolor("#c9c9c9")
        cell.set_linewidth(0.4)
        cell.PAD = 0.065
        if row_index == 0:
            cell.set_facecolor("#eaf3f7")
            cell.set_text_props(weight="bold")
    axis.text(0, 0.10, "Primary denominator: distinct barcode-accession pairs.\n"
              "HSP-weighted values are secondary. Categories describe\n"
              "record titles, not the function of each matched interval.",
              transform=axis.transAxes, fontsize=10, va="top", linespacing=1.5)


def save_figure(figure, output, stem):
    figure.savefig(output / f"{stem}.png", dpi=300, facecolor="white")
    figure.savefig(output / f"{stem}.pdf", facecolor="white")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    arguments = parser.parse_args()
    arguments.out.mkdir(parents=True, exist_ok=True)
    taxa = load_csv(arguments.data / "Fig2A_taxid_groups.csv")
    rankings = load_csv(arguments.data / "Fig2B_barcode_ranking.csv")
    composition = load_csv(arguments.data / "Fig2C_composition.csv")
    report = json.loads((arguments.data / "audit_summary.json").read_text())
    assert report["panels"]["barcode_only"]["barcode_accession_pairs"] == 363
    assert report["panels"]["motif_barcode"]["strict_hsps"] == 0
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "pdf.fonttype": 42, "axes.labelsize": 11})
    figure = plt.figure(figsize=(10.2, 13.9358))
    draw_taxa(figure.add_axes((0.09, 0.62, 0.40, 0.32)), taxa)
    draw_composition(figure.add_axes((0.05, 0.12, 0.45, 0.30)), composition)
    draw_ranking(figure.add_axes((0.61, 0.095, 0.36, 0.845)), rankings)
    save_figure(figure, arguments.out, "Figure2_v29")
    figure = plt.figure(figsize=(6.8, 6.0))
    draw_taxa(figure.add_axes((0.13, 0.28, 0.82, 0.60)), taxa)
    save_figure(figure, arguments.out, "Figure2A_v29")
    figure = plt.figure(figsize=(6.8, 13.9))
    draw_ranking(figure.add_axes((0.14, 0.10, 0.79, 0.84)), rankings)
    save_figure(figure, arguments.out, "Figure2B_v29")
    print("Rendered Figure 2 composite (A/B plus composition table) and separate A/B PNG/PDF panels")


if __name__ == "__main__":
    main()