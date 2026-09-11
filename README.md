# Demultiplexing and Residual Barcode Analysis

Calculation scripts for comparing Dorado and Barbell preprocessing, auditing complete-query barcode-reference matches, and checking detector sensitivity. The main experimental workflow uses ONT SQK-NBD114-96; the reference screen is a focused case study, not a general database-contamination census or cross-platform benchmark.

## Privacy and Repository Contents

This repository contains calculation/pipeline scripts, regression tests and this README only. It contains no sequencing inputs, donor or sample records, read identifiers, per-experiment measurements, result tables, taxonomy responses, manuscript files or reviews. Users supply their own data and barcode definitions locally.

The pipeline can produce sensitive intermediate files, including FASTQ output, read-level SQLite databases, examples of problematic reads and per-experiment reports. Those local outputs are **not approved for publication** merely because the code is public. Keep them outside Git and follow the applicable data-privacy policy. Do not use `git add .` in a directory containing study data.

A separately distributed `supplementary_source_data.zip` can contain the same release scripts and selected highly processed, cross-experiment summaries. That archive is not stored in this repository. A Zenodo version DOI should identify an archived release of this code; no DOI is asserted here before it has been issued.

## Requirements

- Python 3.12 or later, with `regex`, `numpy`, `scipy`, `pandas` and `matplotlib`.
- Dorado 1.1.1 and Dorado 2.0.0 for the two basecalling configurations. The pipeline uses a `dorado200` conda environment for the latter by default; override with `--conda-env-200`.
- Barbell 0.3.2, samtools 1.19.2 and NanoStat 1.6.0 for the recorded preprocessing workflow. Supply executable paths when they are not on PATH.
- ONT POD5 tooling if FAST5 conversion is needed; the pipeline takes POD5 input.
- Optional BLAST+ (`blastn`, `blastdbcmd`) and a locally prepared nucleotide database for barcode screening.
- Optional R with `ggplot2`, `readr` and Cairo graphics support for the current GraphPad-style Figure 2 renderer. It generates PNG/PDF files, not GraphPad project files.
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

In Barbell 0.3.2, `--maximize` selects additional accepted configuration patterns at the filter step after unchanged annotation. It does not increase barcode-search sensitivity or enable `--use-extended`. Retention and residual-pattern removal do not establish correct sample assignment or clinical safety.

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

The BLAST driver requires a user-supplied barcode CSV, output directory and database prefix. If the raw exports already exist, start with the audit command; do not rerun the search unnecessarily.

```bash
THREADS=8 bash scripts/run_barcode_core_nt_blast.sh \
  /path/to/barcodes.csv /path/to/local_blast_output /path/to/core_nt/core_nt

.venv/bin/python scripts/audit_complete_blast.py \
  --barcode-only /path/to/local_blast_output/barcodes_only_blast_raw.tsv \
  --motif-barcode /path/to/local_blast_output/barcodes_with_ont_adapter_blast_raw.tsv \
  --out /path/to/audited_figure2
```

The audit requires `pident = 100`, `mismatch = 0`, `gapopen = 0`, query endpoints spanning 1 to `qlen`, alignment length equal to `qlen`, and an ungapped subject span of the same length. Query length and subject coordinates must be positive. `qcovs` is subject-level coverage and is not proof of individual-HSP completeness. Raw exports use the 17-column schema defined in the audit script.

Distinct barcode-accession pairs are primary; HSPs are secondary. Grouping retains exported TaxIDs rather than parsing names from record titles. Record-title categories do not independently establish the function or biological origin of a matched interval. The audit writes detailed local HSP/accession working tables as well as aggregate plot tables; keep the detailed outputs out of Git.

### Taxonomy and GraphPad-Style Figure 2

Map exported TaxIDs to NCBI Taxonomy scientific names, ranks and lineages, then render:

```bash
.venv/bin/python scripts/map_figure2_taxonomy.py \
  --input /path/to/audited_figure2/Fig2A_taxid_groups.csv \
  --out /path/to/audited_figure2 --cache /path/to/taxonomy_cache

Rscript scripts/plot_audited_figure2_graphpad.R \
  --data /path/to/audited_figure2 --out /path/to/rendered_figure2
```

The first mapping call sends only the TaxIDs to NCBI. Subsequent calls reuse the exact cached XML and request metadata when available, checking its checksum. Missing/duplicate mappings fail explicitly; aliases are recorded without merging the original exported groups. Names and ranks are current NCBI annotations, not independent biological validation or a reconstruction of taxonomy at the BLAST database-build date. Informational species-ancestor columns are not used for regrouping.

The current R renderer preserves the GraphPad-style axes, typography and palette, with named/ranked TaxID labels and connector lines in panel A, accession counts in B, and a composition table in C. It is a study-specific reproduction script: expected totals and selected label positions are deliberately checked. For another dataset, use the general audit/mapping outputs and adapt the figure assertions and layout deliberately. The author's separately assembled Prism stacked panel C is not recreated by this table-C renderer.

The previous Python renderer remains for historical reproduction of the unnamed layout. The superseded title-parsing/risk-class Figure 2 helpers have been retired; use the strict audit and named R renderer for the corrected workflow.

The motif panel prefixes `TTTCTGTTGGTGCTGATATTGCT` separately to the forward barcode and its reverse complement. It is not a complete native ligation construct and its reverse complement. BLAST already searches both strands. The driver's top-ten exports are capped; the strict audit reads all raw exported rows. Differences between query constructions do not isolate adapter context from length or establish suppression of biological cross-reactivity.

## Pooled Read-Length and Yield Tables

Use the integer length-bin exports from all experiments and the saved yield summary:

```bash
.venv/bin/python scripts/make_s1_prism_tables.py \
  --bins /path/to/experiment_a/length_distribution_bins.csv \
         /path/to/experiment_b/length_distribution_bins.csv \
  --yield-summary /path/to/summary_a_read_yield.csv \
  --out /path/to/figureS1_tables
```

Supply every experiment's bin file, not just the two illustrative paths. The exporter reconciles bin totals with yields and writes mean/SD/N yield summaries, pooled counts/percentages and density-format XY tables. Unequal-width bins need density for a continuous histogram; the open-ended overflow bin has no defined width. Aggregate mean/SD/N cannot recreate individual experiment symbols. Logarithmic axis formatting does not change the underlying percentages or their denominators.

## Bounded Detector Validation

The occurrence-enumerating detector is a separate validation tool; it does not replace the original production scanner. It retains the original edit budget, checks candidate spans and resolves overlaps with a deterministic greedy policy. Adjacent intervals are retained; ambiguous positive-length overlaps are rejected.

```bash
.venv/bin/python scripts/validate_detector_occurrences.py \
  --barcode-csv /path/to/barcodes.csv \
  --experiments /path/to/experiment_a_demux /path/to/experiment_b_demux \
  --samples-per-condition 64 --windows 40 80 120 --seed 20260911 \
  --synthetic --out /path/to/private_validation

.venv/bin/python scripts/summarize_detector_validation.py \
  --source /path/to/private_validation --out /path/to/aggregate_validation
```

Supply all intended experiment directories. Each must contain the condition output folders expected by the CLI; include only actual output FASTQs, not duplicated input/merged intermediate files. Sampling is bounded, restricted to reads shorter than 5,000 bp and paired by experiment, not read ID. The random-byte proposal uses approximate inverse-record-span weighting. Synthetic-only validation is available by omitting `--experiments`.

The aggregate exporter separates study-pilot, experiment-paired, synthetic-recovery and synthetic-background summaries. These are not replacements for full-study length distributions or residual-rate estimates. Model backgrounds and mononucleotide shuffles are not matched biological barcode-free controls; zero observed events do not demonstrate zero background, equal sensitivity or rare-event precision. Files labelled `_LOCAL` and all sequences/identifiers must remain private.

## Tests

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python -m unittest discover -s scripts/tests -v
```

The 20 tests cover strict HSP criteria, taxonomy parsing/alias and failure handling, and occurrence-detector behaviour including comparison with the production best-hit arm. They do not require study reads, a BLAST database or network taxonomy access.

## Measurement Limits and Release Notes

- The regex budget is at most one substitution, one insertion and one deletion per barcode, not any three edits. One best match is retained per barcode orientation. Repeated copies are not exhaustively enumerated; a terminal best match can hide an internal copy, and one short-read match can satisfy both terminal-window tests.
- Condition-level record ratios, barcode-ID multiplicity and location summaries are detector-defined endpoints, not calibrated sample-crosstalk or clinical false-positive rates.
- The strict per-HSP audit supersedes the historical `qcovs`-based filter. The legacy BLAST driver's convenience summaries are not substitutes for the audited outputs.
- The legacy `write_stat_tests` branch in the scanner queries `dorado`/`barbell` labels. Its output is not a valid pairwise-test table for the c0-c6 design. Use the separate paired experiment-level analysis; the legacy branch was not silently redefined for the published comparison.
- `stat_mean` length-bin summaries are unweighted experiment means. Pooled percentages require weighting by each condition's output read counts. Geometric summaries are not additive.
- Barcode-like sequences in experimental reads and matching kit barcode sequences in references are complementary observations. They do not establish that the residual reads caused false taxonomic assignments or generated the surveyed reference records.
- NCBI name/rank annotation does not turn mixed-rank bacterial, eukaryotic and viral groups into a bacteria-only species census; title-keyword categories do not establish interval function.

This update publishes the portable scripts and tests from the current supplementary code payload, including the strict audit, pooled figure tables, detector pilot, verified taxonomy mapping and connector-line GraphPad renderer. The existing production scanner remains unchanged. No sequencing, full-study rescan, new BLAST search or clinical classification was run during publication preparation. The tested environment is not asserted to reproduce the exact historical production environment.

For a reproducible citation, record the Git commit used and cite the Zenodo **version DOI** after the matching release has been archived. A concept DOI does not pin a specific release.