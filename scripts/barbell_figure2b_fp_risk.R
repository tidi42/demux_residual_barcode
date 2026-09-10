#!/usr/bin/env Rscript
###############################################################################
# barbell_figure2b_fp_risk.R
#
# Figure 2B - Barcode vulnerability ranking, GraphPad-Prism style. Same graph
# type as the original: horizontal, sorted bars (barcode-only condition),
# coloured by risk class and annotated with the number of distinct organisms
# ("N spp."). One figure is produced.
#
# Strict perfect-hit filter throughout: mismatch = 0 AND pident = 100 AND
# qcovs = 100 (subject-level coverage, not per-HSP completeness).
#
# 2026-09-08: switched from the ten-best-hits-per-orientation export to the
# COMPLETE BLAST output. The capped export imposed a ceiling of 10 x 2 = 20
# hits per barcode, at which sixteen barcodes sat tied; those flat bars were an
# artefact of the export, not a measurement. Bars now span 2-218 hits.
#
# Input  (figures/fig2)    : Fig2B_barcode_ranking_FULLDATA.csv
#                            (validated against barcodes_only_blast_raw.tsv:
#                             1,252 perfect HSPs / 366 distinct barcode-subject
#                             pairs / 289 organisms over 55 barcodes)
# Output (figures/fig2)    : Fig2B_barcode_vulnerability_ranking_graphpad.png
#                            Fig2B_barcode_vulnerability_ranking_graphpad.csv
###############################################################################

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/barbell_figure2b_fp_risk.R [TABLE_DIRECTORY]\n")
  quit(status = 0)
}
suppressPackageStartupMessages(library(tidyverse))
if (capabilities("cairo")) options(bitmapType = "cairo")

# --- paths (relative to project root) ---------------------------------------
outdir   <- if (length(arguments)) arguments[1] else "figures/fig2"
full_csv <- file.path(outdir, "Fig2B_barcode_ranking_FULLDATA.csv")
png_out  <- file.path(outdir, "Fig2B_barcode_vulnerability_ranking_graphpad.png")
csv_out  <- file.path(outdir, "Fig2B_barcode_vulnerability_ranking_graphpad.csv")

# --- colour-blind-safe warm palette (Okabe-Ito): yellow -> orange -> red ----
RISK_COLS <- c("Low (1-9)"        = "#F0E442",   # yellow
               "Moderate (10-19)" = "#E69F00",   # orange
               "High (>=20 hits)" = "#D55E00")   # vermillion (red-orange)

# --- GraphPad-Prism-like theme ----------------------------------------------
theme_prism_like <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line          = element_line(colour = "black", linewidth = 0.9, lineend = "square"),
      axis.ticks         = element_line(colour = "black", linewidth = 0.9),
      axis.ticks.length  = unit(0.18, "cm"),
      axis.text          = element_text(colour = "black"),
      axis.title         = element_text(colour = "black", face = "bold"),
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.title.position = "plot",
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      legend.key         = element_blank(),
      panel.grid         = element_blank()
    )
}

# --- load the complete-output per-barcode table -----------------------------
plot_df <- read_csv(full_csv, show_col_types = FALSE) %>%
  filter(panel == "Barcode only") %>%
  transmute(barcode,
            n_perfect    = n_perfect_hsp,
            n_subjects   = n_perfect_distinct_subjects,
            n_organisms,
            n_species) %>%
  arrange(desc(n_perfect)) %>%
  mutate(
    barcode    = fct_reorder(barcode, n_perfect),
    risk_class = factor(case_when(
      n_perfect >= 20 ~ "High (>=20 hits)",
      n_perfect >= 10 ~ "Moderate (10-19)",
      TRUE            ~ "Low (1-9)"),
      levels = names(RISK_COLS))
  )

stopifnot(nrow(plot_df) > 0, all(plot_df$n_perfect >= plot_df$n_subjects))
write_csv(arrange(plot_df, desc(n_perfect)), csv_out)

# --- plot: sorted horizontal bars + "N spp." labels -------------------------
p <- ggplot(plot_df, aes(n_perfect, barcode, fill = risk_class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(n_species, " spp.")),
            hjust = -0.15, size = 6.15, colour = "grey20") +
  scale_fill_manual(values = RISK_COLS, name = "Risk class", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Perfect full-query BLAST hits", y = NULL) +
  theme_prism_like(base_size = 18.7) +
  theme(axis.line    = element_line(colour = "black", linewidth = 1.8, lineend = "square"),
        axis.ticks   = element_line(colour = "black", linewidth = 1.2),
        axis.text.y  = element_text(size = 17.6, face = "bold"),
        axis.text.x  = element_text(size = 18.7, face = "bold"),
        axis.title.x = element_text(size = 20.9),
        legend.title = element_text(size = 21.6, face = "bold"),
        legend.text  = element_text(size = 20.4))

n_bc <- nrow(plot_df)
ggsave(png_out, p, width = 9.5, height = max(5, n_bc * 0.501 + 3), dpi = 600)
message("Saved: ", png_out, "  (", n_bc, " barcodes, ",
        min(plot_df$n_perfect), "-", max(plot_df$n_perfect), " hits)")
message("Saved: ", csv_out)
