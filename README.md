# Demultiplexing and Residual Barcode Analysis

Calculation scripts for comparing Dorado and Barbell preprocessing and for an exploratory barcode-to-database similarity screen.

## Privacy and Repository Contents

This repository contains calculation/pipeline scripts and this README only. It contains no sequencing inputs, donor or sample records, read identifiers, per-experiment measurements, analysis results, manuscript files or reviews. Users supply their own data and barcode definitions locally.

The pipeline can produce sensitive intermediate files, including FASTQ output, read-level SQLite databases, examples of problematic reads and per-experiment reports. Those local outputs are **not approved for publication** merely because the code is public. Keep them outside Git and follow the applicable data-privacy policy. Do not use `git add .` in a directory containing study data.

A separately distributed `supplementary_source_data.zip` can contain the same release scripts and selected highly processed, cross-experiment summaries. That archive is not stored in this repository. A Zenodo version DOI should identify an archived release of this code; no DOI is asserted here before it has been issued.

## Requirements

- Python 3.12 or later, with `regex`, `numpy`, `scipy`, `pandas` and `matplotlib`.
- Dorado 1.1.1 and Dorado 2.0.0 for the two basecalling configurations. The pipeline uses a `dorado200` conda environment for the latter by default; override with `--conda-env-200`.
- Barbell 0.3.2, samtools 1.19.2 and NanoStat 1.6.0 for the recorded preprocessing workflow. Supply executable paths when they are not on PATH.
- ONT POD5 tooling if FAST5 conversion is needed; the pipeline takes POD5 input.
- Optional BLAST+ (`blastn`, `blastdbcmd`) and a locally prepared nucleotide database for barcode screening.
- Optional R with `tidyverse` for the two database-summary plotting scripts. These are plotting helpers, not GraphPad project files.
- Optional minimap2 only for the legacy auxiliary barcode scan; it is not the approximate-regex endpoint described below.

Use an isolated environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install regex numpy scipy pandas matplotlib
```

The supported CLI entry points provide `--help`. Test the external tools independently before a full run. Thread counts and GPU requirements should match the host; the recorded 120-thread setting is not a recommended minimum.

## Preprocessing Conditions

| Condition | Configuration |
| --- | --- |
| C0 | Dorado 1.1.1, HAC v5.0.0, no demultiplexing or trimming |
| C1 | Dorado 2.0.0, HAC v6.0.0, no demultiplexing or trimming |
| C2/C3 | Corresponding Dorado kit demultiplexing and trimming |
| C4/C5 | Barbell kit mode with `--maximize`, applied to C0/C1 |
| C6 | Barbell applied to the concatenated C3 output; a selected diagnostic subset |

The release and basecalling model change together. C4/C5 are parallel alternatives to C2/C3, not sequential processing of them. C6 is not a matched whole-population comparator.

## Pipeline

Run from the repository root. The barcode CSV is supplied by the user and must have columns `Component`, `Forward sequence` and `Reverse sequence`. No barcode-input file is bundled.

```bash
.venv/bin/python scripts/barbell_paper.py \
  --pod5-dir /path/to/retained/pod5 \
  --out outputs/experiment_a_demux \
  --barcode-csv /path/to/barcodes.csv \
  --dorado-111 /path/to/dorado-1.1.1/bin/dorado \
  --conda-env-200 dorado200 \
  --barbell /path/to/barbell \
  --samtools /path/to/samtools \
  --python .venv/bin/python \
  --threads 8
```

Paths in this release are portable defaults or explicit arguments; they do not encode the source study's experiment labels. Input selection, including adaptive-sampling pass/fail handling, is performed before invoking the pipeline and must be appropriate to the user's data.

For already prepared FASTQs, the scanner can be invoked directly:

```bash
.venv/bin/python scripts/compare_dorado_barbell_outputs.py \
  --condition-dirs /path/to/c0 /path/to/c1 /path/to/c2 /path/to/c3 /path/to/c4 /path/to/c5 /path/to/c6 \
  --condition-names c0 c1 c2 c3 c4 c5 c6 \
  --barcode-csv /path/to/barcodes.csv \
  --out outputs/experiment_a_demux/cmp_all \
  --nanostat /path/to/NanoStat \
  --skip-barcode-scan \
  --threads 8
```

Here `--skip-barcode-scan` skips the legacy **minimap2** auxiliary scan, not the regex endpoint. The regex scan is enabled by `--barcode-csv`; omitting that CSV skips residual-pattern detection. Choose a new output directory for a new analysis. `--overwrite` replaces an existing comparison database and should be used only deliberately.

## Aggregation and Paired Statistics

The local directory pattern is `outputs/{experiment}_demux/cmp_all/`. Supply the experiment names explicitly; the examples below are placeholders, not subject identifiers.

```bash
.venv/bin/python scripts/barbell_figures_summa.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures --skip-figures

.venv/bin/python scripts/barbell_figure_f_paired.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures/figure_f_violin

.venv/bin/python scripts/make_barcode_multiplicity_csvs.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures/graphpad
```

`barbell_figures.py` generates per-condition summaries from one comparison database. `barbell_figure_f_violin.py` calculates location-rate summaries and optional rarefaction displays from supplied experiment directories. `make_summary_f_log.py` takes `--input`, `--experiments` and `--out` to generate log-scale location summaries. Their generated replicate-level files are local analysis outputs, not privacy-filtered publication data.

The primary inferential endpoint is the paired, experiment-level internal-match rate in reads per million output reads. The paired script reports unadjusted Wilcoxon signed-rank tests and ratios. Its legacy plotting defaults are not a claim about the formatting of any final manuscript figure.

## Barcode-Database Screen

The BLAST driver requires a user-supplied barcode CSV, output directory and database prefix:

```bash
THREADS=8 bash scripts/run_barcode_core_nt_blast.sh \
  /path/to/barcodes.csv /path/to/local_blast_output /path/to/core_nt/core_nt

.venv/bin/python scripts/make_fig2_graphpad_tables.py \
  --barcode-hits /path/to/local_blast_output/barcodes_only_blast_raw.tsv \
  --motif-hits /path/to/local_blast_output/barcodes_with_ont_adapter_blast_raw.tsv \
  --out figures/fig2

.venv/bin/python scripts/keyword_precedence_sensitivity.py \
  --barcode-hits /path/to/local_blast_output/barcodes_only_blast_raw.tsv \
  --motif-hits /path/to/local_blast_output/barcodes_with_ont_adapter_blast_raw.tsv
```

The table generator writes complete HSP-count, subject/title-group and record-annotation tables. Optional R helpers take the table directory as their first argument:

```bash
Rscript scripts/barbell_figure2a_shared_vulnerability_bubble.R figures/fig2
Rscript scripts/barbell_figure2b_fp_risk.R figures/fig2
```

The motif panel prefixes `TTTCTGTTGGTGCTGATATTGCT` separately to the forward barcode and its reverse complement. It is not a complete native ligation construct and its reverse complement. BLAST already searches both strands. The driver's top-ten exports are capped; the complete-table generator instead reads the raw BLAST output.

## Measurement Limits and Release Notes

- The regex budget is at most one substitution, one insertion and one deletion per barcode, not any three edits. One best match is retained per barcode orientation. Repeated copies are not exhaustively enumerated; a terminal best match can hide an internal copy, and one short-read match can satisfy both terminal-window tests.
- Condition-level record ratios, barcode-ID multiplicity and location summaries are detector-defined endpoints, not calibrated sample-crosstalk or clinical false-positive rates.
- The historical BLAST filter is `mismatch = 0`, `pident = 100`, `qcovs = 100`. `qcovs` is query coverage per subject across HSPs, not the span of an individual HSP. No full-query validation or new alignment totals are implied by this release.
- The legacy `write_stat_tests` branch in the scanner queries `dorado`/`barbell` labels. Its output is not a valid pairwise-test table for the c0-c6 design. Use the separate paired experiment-level analysis; the legacy branch was not silently redefined for the published comparison.
- `stat_mean` length-bin summaries are unweighted experiment means. Pooled percentages require weighting by each condition's output read counts. Geometric summaries are not additive.
- Title parsing and keyword precedence yield descriptive database-record groups/categories, not independent taxonomic or functional validation.

This release makes executable paths and experiment lists configurable, builds the complete database-summary tables from explicit inputs, and removes fixed study-total assertions. The regex detector and operational BLAST-filter rules are unchanged. Python/R/shell syntax and CLI entry points were checked; synthetic tests covered detector equivalence, paired statistics, HSP filtering/de-duplication and percentage-table totals. Full basecalling, BLAST searches and private-study recomputation were not run during release preparation. The tested environment is not asserted to be an export of the historical production environment.

For a reproducible citation, record the Git commit used and cite the Zenodo **version DOI** after the matching release has been archived. A concept DOI does not pin a specific release.