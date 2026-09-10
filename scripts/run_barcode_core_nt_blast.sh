#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Barcode BLAST analysis against core_nt
#
# Input:
#   User-supplied barcode CSV (first argument)
#
# Output:
#   User-selected local output directory (second argument)
#
# Requirements:
#   blastn installed
#   core_nt BLAST database available
#   python3 available
###############################################################################

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    printf '%s\n' 'Usage: bash scripts/run_barcode_core_nt_blast.sh BARCODE_CSV OUTPUT_DIR BLAST_DATABASE' 'Optional environment: THREADS (default 120)'
    exit 0
fi
INPUT_CSV="${1:?Supply BARCODE_CSV, OUTPUT_DIR and BLAST_DATABASE}"
OUTDIR="${2:?Supply OUTPUT_DIR}"
DB="${3:?Supply the BLAST database prefix}"
THREADS="${THREADS:-120}"

# ONT adapter sequence used for adapter-plus-barcode search.
# This is the common ONT adapter motif often seen in ligation/native workflows.
ONT_ADAPTER="TTTCTGTTGGTGCTGATATTGCT"

# Number of raw BLAST hits collected before post-ranking.
# We collect more than 10 because many short barcode hits may have 100% identity.
MAX_RAW_HITS=5000

# Final number of hits retained per query sequence.
TOP_N=10

mkdir -p "$OUTDIR"

BARCODE_FASTA="$OUTDIR/barcodes_only.fasta"
ADAPTER_FASTA="$OUTDIR/barcodes_with_ont_adapter.fasta"

BARCODE_BLAST_RAW="$OUTDIR/barcodes_only_blast_raw.tsv"
ADAPTER_BLAST_RAW="$OUTDIR/barcodes_with_ont_adapter_blast_raw.tsv"

BARCODE_TOP10="$OUTDIR/barcodes_only_top10_hits.csv"
ADAPTER_TOP10="$OUTDIR/barcodes_with_ont_adapter_top10_hits.csv"

BARCODE_OVERVIEW="$OUTDIR/barcodes_only_analysis_overview.csv"
ADAPTER_OVERVIEW="$OUTDIR/barcodes_with_ont_adapter_analysis_overview.csv"

COMBINED_OVERVIEW="$OUTDIR/combined_analysis_overview.csv"

echo "[INFO] Input CSV: $INPUT_CSV"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] BLAST database: $DB"

###############################################################################
# 1) Generate FASTA query files from CSV
###############################################################################

python3 - <<PY
import csv
from pathlib import Path

input_csv = Path("$INPUT_CSV")
barcode_fasta = Path("$BARCODE_FASTA")
adapter_fasta = Path("$ADAPTER_FASTA")
adapter = "$ONT_ADAPTER"

def clean_seq(seq):
    return seq.strip().upper().replace(" ", "").replace("\t", "")

with input_csv.open(newline="") as f, \
     barcode_fasta.open("w") as out1, \
     adapter_fasta.open("w") as out2:

    reader = csv.DictReader(f)

    required = {"Component", "Forward sequence", "Reverse sequence"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"Missing required CSV columns: {missing}")

    for row in reader:
        component = row["Component"].strip()
        fwd = clean_seq(row["Forward sequence"])
        rev = clean_seq(row["Reverse sequence"])

        if not component or not fwd or not rev:
            continue

        # Barcode-only queries
        out1.write(f">{component}_forward\\n{fwd}\\n")
        out1.write(f">{component}_reverse\\n{rev}\\n")

        # Adapter-plus-barcode queries
        # Prefix mode: models adapter directly upstream of barcode.
        out2.write(f">{component}_forward_adapter_plus_barcode\\n{adapter}{fwd}\\n")
        out2.write(f">{component}_reverse_adapter_plus_barcode\\n{adapter}{rev}\\n")

print(f"Created: {barcode_fasta}")
print(f"Created: {adapter_fasta}")
PY

###############################################################################
# 2) Check files and database access
###############################################################################

if [[ ! -s "$BARCODE_FASTA" ]]; then
    echo "[ERROR] Barcode FASTA is empty: $BARCODE_FASTA" >&2
    exit 1
fi

if [[ ! -s "$ADAPTER_FASTA" ]]; then
    echo "[ERROR] Adapter FASTA is empty: $ADAPTER_FASTA" >&2
    exit 1
fi

if ! command -v blastn >/dev/null 2>&1; then
    echo "[ERROR] blastn not found in PATH." >&2
    exit 1
fi

echo "[INFO] Checking BLAST database..."
blastdbcmd -db "$DB" -info >/dev/null

###############################################################################
# 3) Run BLAST
###############################################################################

# Important:
# -task blastn-short is suitable for short barcode-like queries.
# -max_target_seqs is intentionally high.
# Final top 10 are selected later using deterministic sorting.
# qcovs is coverage aggregated per subject, not individual-HSP completeness.

BLAST_OUTFMT="6 qseqid saccver staxids sscinames pident length qlen qcovs mismatch gapopen qstart qend sstart send evalue bitscore stitle"

echo "[INFO] Running BLAST for barcode-only queries..."

blastn \
    -task blastn-short \
    -query "$BARCODE_FASTA" \
    -db "$DB" \
    -out "$BARCODE_BLAST_RAW" \
    -outfmt "$BLAST_OUTFMT" \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_RAW_HITS" \
    -dust no \
    -soft_masking false

echo "[INFO] Running BLAST for adapter-plus-barcode queries..."

blastn \
    -task blastn-short \
    -query "$ADAPTER_FASTA" \
    -db "$DB" \
    -out "$ADAPTER_BLAST_RAW" \
    -outfmt "$BLAST_OUTFMT" \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_RAW_HITS" \
    -dust no \
    -soft_masking false

###############################################################################
# 4) Rank hits and create overviews
###############################################################################

python3 - <<PY
import csv
from pathlib import Path
from collections import defaultdict

top_n = int("$TOP_N")

columns = [
    "qseqid",
    "saccver",
    "staxids",
    "sscinames",
    "pident",
    "length",
    "qlen",
    "qcovs",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
]

numeric_float = {"pident", "qcovs", "evalue", "bitscore"}
numeric_int = {"length", "qlen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"}

def read_fasta_ids(path):
    ids = []
    with Path(path).open() as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def parse_blast_tsv(path):
    hits = []
    path = Path(path)

    if not path.exists() or path.stat().st_size == 0:
        return hits

    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < len(columns):
                continue

            rec = dict(zip(columns, row[:len(columns)]))

            for key in numeric_float:
                try:
                    rec[key] = float(rec[key])
                except ValueError:
                    rec[key] = None

            for key in numeric_int:
                try:
                    rec[key] = int(rec[key])
                except ValueError:
                    rec[key] = None

            hits.append(rec)

    return hits

def sort_key(hit):
    # Deterministic ranking:
    # 1. highest percent identity
    # 2. highest query coverage
    # 3. longest alignment
    # 4. highest bit score
    # 5. lowest e-value
    # 6. stable accession/title ordering
    return (
        -(hit["pident"] if hit["pident"] is not None else -1),
        -(hit["qcovs"] if hit["qcovs"] is not None else -1),
        -(hit["length"] if hit["length"] is not None else -1),
        -(hit["bitscore"] if hit["bitscore"] is not None else -1),
        hit["evalue"] if hit["evalue"] is not None else float("inf"),
        hit["saccver"],
        hit["stitle"],
    )

def write_top10_and_overview(raw_tsv, fasta, top10_csv, overview_csv, analysis_type):
    query_ids = read_fasta_ids(fasta)
    hits = parse_blast_tsv(raw_tsv)

    by_query = defaultdict(list)
    for h in hits:
        by_query[h["qseqid"]].append(h)

    top_rows = []
    overview_rows = []

    for qid in query_ids:
        qhits = by_query.get(qid, [])
        # Keep only hits with at most 2 mismatches
        qhits = [h for h in qhits if h.get("mismatch") is not None and h["mismatch"] <= 2]
        qhits_sorted = sorted(qhits, key=sort_key)
        selected = qhits_sorted[:top_n]

        for rank, h in enumerate(selected, start=1):
            row = {
                "analysis_type": analysis_type,
                "rank": rank,
                **h,
            }
            top_rows.append(row)

        if selected:
            best = selected[0]
            hits_100 = sum(
                1 for h in qhits
                if h["pident"] == 100.0
            )
            hits_100_full_query = sum(
                1 for h in qhits
                if h["pident"] == 100.0 and h["qcovs"] == 100.0
            )

            overview_rows.append({
                "analysis_type": analysis_type,
                "qseqid": qid,
                "number_of_raw_hits": len(qhits),
                "number_of_100_percent_identity_hits": hits_100,
                "number_of_100_percent_identity_and_100_percent_qcov_hits": hits_100_full_query,
                "best_saccver": best["saccver"],
                "best_staxids": best["staxids"],
                "best_sscinames": best["sscinames"],
                "best_pident": best["pident"],
                "best_qcovs": best["qcovs"],
                "best_alignment_length": best["length"],
                "best_evalue": best["evalue"],
                "best_bitscore": best["bitscore"],
                "best_stitle": best["stitle"],
            })
        else:
            overview_rows.append({
                "analysis_type": analysis_type,
                "qseqid": qid,
                "number_of_raw_hits": 0,
                "number_of_100_percent_identity_hits": 0,
                "number_of_100_percent_identity_and_100_percent_qcov_hits": 0,
                "best_saccver": "",
                "best_staxids": "",
                "best_sscinames": "",
                "best_pident": "",
                "best_qcovs": "",
                "best_alignment_length": "",
                "best_evalue": "",
                "best_bitscore": "",
                "best_stitle": "",
            })

    top_fields = ["analysis_type", "rank"] + columns
    overview_fields = [
        "analysis_type",
        "qseqid",
        "number_of_raw_hits",
        "number_of_100_percent_identity_hits",
        "number_of_100_percent_identity_and_100_percent_qcov_hits",
        "best_saccver",
        "best_staxids",
        "best_sscinames",
        "best_pident",
        "best_qcovs",
        "best_alignment_length",
        "best_evalue",
        "best_bitscore",
        "best_stitle",
    ]

    with Path(top10_csv).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=top_fields)
        writer.writeheader()
        writer.writerows(top_rows)

    with Path(overview_csv).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=overview_fields)
        writer.writeheader()
        writer.writerows(overview_rows)

    return overview_rows

barcode_overview = write_top10_and_overview(
    "$BARCODE_BLAST_RAW",
    "$BARCODE_FASTA",
    "$BARCODE_TOP10",
    "$BARCODE_OVERVIEW",
    "barcode_only"
)

adapter_overview = write_top10_and_overview(
    "$ADAPTER_BLAST_RAW",
    "$ADAPTER_FASTA",
    "$ADAPTER_TOP10",
    "$ADAPTER_OVERVIEW",
    "adapter_plus_barcode"
)

combined = barcode_overview + adapter_overview

overview_fields = [
    "analysis_type",
    "qseqid",
    "number_of_raw_hits",
    "number_of_100_percent_identity_hits",
    "number_of_100_percent_identity_and_100_percent_qcov_hits",
    "best_saccver",
    "best_staxids",
    "best_sscinames",
    "best_pident",
    "best_qcovs",
    "best_alignment_length",
    "best_evalue",
    "best_bitscore",
    "best_stitle",
]

with Path("$COMBINED_OVERVIEW").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=overview_fields)
    writer.writeheader()
    writer.writerows(combined)

print("Created top10 and overview files.")
PY

###############################################################################
# 5) Done
###############################################################################

echo ""
echo "[DONE] Barcode BLAST analysis finished."
echo ""
echo "Output files:"
echo "  $BARCODE_FASTA"
echo "  $ADAPTER_FASTA"
echo "  $BARCODE_BLAST_RAW"
echo "  $ADAPTER_BLAST_RAW"
echo "  $BARCODE_TOP10"
echo "  $ADAPTER_TOP10"
echo "  $BARCODE_OVERVIEW"
echo "  $ADAPTER_OVERVIEW"
echo "  $COMBINED_OVERVIEW"
echo ""
