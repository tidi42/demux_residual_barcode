#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
library(grid)

arguments <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag) {
  position <- match(flag, arguments)
  if (is.na(position) || position == length(arguments)) stop("Required argument: ", flag)
  arguments[[position + 1L]]
}
data_dir <- argument_value("--data")
out_dir <- argument_value("--out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_colours <- c("Barcode only" = "#D55E00", "Motif + barcode" = "#0072B2")
theme_prism_like <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.9, lineend = "square"),
      axis.ticks = element_line(colour = "black", linewidth = 0.9),
      axis.ticks.length = unit(0.18, "cm"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.title.position = "plot",
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid = element_blank()
    )
}

taxa <- read_csv(file.path(data_dir, "Fig2A_taxid_groups.csv"), show_col_types = FALSE)
taxa <- taxa[taxa$panel == "barcode_only", ]
taxa <- taxa[order(-taxa$barcode_accession_pairs, taxa$taxon_group), ]
taxa$condition <- factor("Barcode only", levels = names(condition_colours))
named_taxa <- read_csv(file.path(data_dir, "Fig2A_taxid_groups_named.csv"), show_col_types = FALSE)
matched <- match(taxa$taxon_group, named_taxa$taxon_group)
stopifnot(!anyDuplicated(named_taxa$taxon_group), !anyNA(matched), nrow(named_taxa) == nrow(taxa))
for (column in c("distinct_barcodes", "distinct_accessions", "barcode_accession_pairs", "n_hsp")) {
  stopifnot(identical(taxa[[column]], named_taxa[[column]][matched]))
}
for (column in c("scientific_name", "taxonomic_rank", "domain", "exported_taxid", "species_name")) {
  taxa[[column]] <- named_taxa[[column]][matched]
}
stopifnot(!anyNA(taxa$scientific_name), !anyNA(taxa$taxonomic_rank))
ranking <- read_csv(file.path(data_dir, "Fig2B_barcode_ranking.csv"), show_col_types = FALSE)
stopifnot(sum(ranking$panel == "barcode_only") == 96,
          sum(ranking$panel == "motif_barcode") == 96,
          all(ranking$distinct_accessions[ranking$panel == "motif_barcode"] == 0))
ranking <- ranking[ranking$panel == "barcode_only" & ranking$distinct_accessions > 0, ]
ranking <- ranking[order(-ranking$distinct_accessions, ranking$barcode), ]
ranking$barcode <- factor(ranking$barcode, levels = rev(ranking$barcode))
composition <- read_csv(file.path(data_dir, "Fig2C_composition.csv"), show_col_types = FALSE)
composition <- composition[composition$panel == "barcode_only", ]
pairs <- composition[composition$unit == "barcode_accession_pairs", ]
hsps <- composition[composition$unit == "HSPs", ]
stopifnot(nrow(taxa) == 143, nrow(ranking) == 55,
          sum(ranking$distinct_accessions) == 363,
          sum(taxa$barcode_accession_pairs) == 363,
          sum(pairs$count) == 363, sum(hsps$count) == 1128)

labels <- head(taxa, 4)
label_expression <- function(name, species, rank, taxid) {
  quoted <- function(value) encodeString(value, quote = '"')
  if (!is.na(species) && startsWith(name, species) && !startsWith(name, "uncultured")) {
    title <- paste0("italic(", quoted(species), ")")
    suffix <- trimws(substring(name, nchar(species) + 1))
    if (nzchar(suffix)) title <- paste(title, quoted(suffix), sep = " ~ ")
  } else {
    lines <- strwrap(name, width = 25)
    title <- quoted(lines[[1]])
    if (length(lines) > 1) {
      for (line in lines[-1]) title <- paste0("atop(", title, ",", quoted(line), ")")
    }
  }
  paste0("atop(", title, ",plain(", quoted(paste0(rank, "; TaxID ", taxid)), "))")
}
labels$label <- mapply(label_expression, labels$scientific_name, labels$species_name,
                      labels$taxonomic_rank, labels$exported_taxid)
stopifnot(setequal(labels$taxon_group, c("taxid:5823", "taxid:2301481", "taxid:562", "taxid:573")))
label_positions <- data.frame(taxon_group = c("taxid:5823", "taxid:2301481", "taxid:562", "taxid:573"),
                              label_x = c(1.3, 4.5, 13.5, 10.5),
                              label_y = c(88, 47, 7.5, 28),
                              connector_x = c(1.72, 4.5, 13.5, 10.5),
                              connector_y = c(75, 36, 9.6, 22),
                              label_hjust = c(0, 0, 0.5, 0.5))
labels <- merge(labels, label_positions, by = "taxon_group", sort = FALSE)

panel_a <- ggplot(taxa, aes(distinct_barcodes, distinct_accessions)) +
  geom_segment(data = labels,
               aes(x = distinct_barcodes, y = distinct_accessions,
                   xend = connector_x, yend = connector_y),
               inherit.aes = FALSE, colour = "grey35", linewidth = 0.45,
               lineend = "round", show.legend = FALSE) +
  geom_point(aes(size = barcode_accession_pairs, fill = condition),
             shape = 21, alpha = 0.80, colour = "grey15", stroke = 0.5, show.legend = TRUE) +
  geom_text(data = labels, aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
            size = 3.75, colour = "grey20", vjust = 0.5, parse = TRUE) +
  scale_fill_manual(values = condition_colours, name = NULL, drop = FALSE) +
  scale_size_area(max_size = 12, breaks = c(5, 25, 50, 75),
                  name = "Barcode-accession pairs") +
  scale_x_continuous(breaks = seq(1, 14), limits = c(0.4, 15.5), expand = c(0, 0)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                labels = c("1", "2", "5", "10", "20", "50", "100"),
                limits = c(0.75, 110), expand = c(0, 0)) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 5, alpha = 1)),
         size = guide_legend(order = 2, nrow = 1)) +
  labs(x = "Distinct barcode IDs per TaxID group",
       y = "Distinct accessions per TaxID group",
      caption = "143 exported groups; names/ranks: NCBI Taxonomy, 11 Sep 2026.\nCoincident points overlap. Motif + barcode: no complete matches.") +
  theme_prism_like(base_size = 16.5) +
  theme(legend.box = "vertical", legend.spacing.y = unit(2, "pt"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        plot.caption = element_text(size = 12, hjust = 0),
        plot.margin = margin(6, 12, 6, 6))

panel_b <- ggplot(ranking, aes(distinct_accessions, barcode)) +
  geom_col(width = 0.7, fill = condition_colours[["Barcode only"]]) +
  geom_text(aes(label = distinct_accessions), hjust = -0.15,
            size = 6.15, colour = "grey20") +
  scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 82), expand = c(0, 0)) +
  labs(x = "Distinct subject accessions per barcode", y = NULL,
       caption = "55/96 barcodes have complete matches; 41 have zero.\nMotif + barcode: zero for all 96 barcodes.") +
  theme_prism_like(base_size = 18.7) +
  theme(axis.line = element_line(colour = "black", linewidth = 1.8, lineend = "square"),
        axis.ticks = element_line(colour = "black", linewidth = 1.2),
        axis.text.y = element_text(size = 17.6, face = "bold"),
        axis.text.x = element_text(size = 18.7, face = "bold"),
        axis.title.x = element_text(size = 20.9),
        plot.caption = element_text(size = 14, hjust = 0),
        plot.margin = margin(10, 12, 10, 6))

save_panel <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
      width = width, height = height, dpi = 600, bg = "white",
      device = grDevices::png, type = "cairo")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
}
save_panel(panel_a, "Figure2A_v29", 9.5, 7)
save_panel(panel_b, "Figure2B_v29", 9.5, max(5, nrow(ranking) * 0.501 + 3))

draw_composition <- function() {
  grid.text("Record-annotation composition", x = 0, y = 0.96,
            just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold"))
  columns <- c(0.01, 0.52, 0.79)
  grid.text(c("Title category", "Pairs\n(n = 363)", "HSPs\n(n = 1,128)"),
            x = columns, y = 0.84, just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.lines(x = c(0, 1), y = c(0.78, 0.78), gp = gpar(lwd = 1.5))
  for (row_index in seq_len(nrow(pairs))) {
    row <- pairs[row_index, ]
    hsp <- hsps[hsps$category == row$category, ]
    row_y <- 0.71 - (row_index - 1) * 0.10
    grid.text(c(sub("Chromosome/genome", "Chromosome/\ngenome", row$category, fixed = TRUE),
                sprintf("%d (%.1f%%)", row$count, row$percent),
                sprintf("%d (%.1f%%)", hsp$count, hsp$percent)),
              x = columns, y = row_y, just = "left", gp = gpar(fontsize = 8.5))
  }
  grid.lines(x = c(0, 1), y = c(0.15, 0.15), gp = gpar(lwd = 1.5))
  grid.text("Primary denominator: barcode-accession pairs.\nCategories describe record titles, not interval function.",
            x = 0, y = 0.10, just = c("left", "top"), gp = gpar(fontsize = 8))
}

draw_composite <- function() {
  grid.newpage()
  compact_a <- panel_a + theme(
    text = element_text(size = 10), axis.text = element_text(size = 9),
    axis.title = element_text(size = 10), legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    plot.caption = element_text(size = 8), plot.margin = margin(4, 4, 4, 4))
  compact_b <- panel_b + theme(
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 9), plot.caption = element_text(size = 8))
  compact_a$layers[[3]]$aes_params$size <- 2.6
  compact_b$layers[[2]]$aes_params$size <- 3
  print(compact_a, vp = viewport(x = 0.265, y = 0.76, width = 0.51, height = 0.40))
  print(compact_b, vp = viewport(x = 0.755, y = 0.50, width = 0.47, height = 0.92))
  pushViewport(viewport(x = 0.27, y = 0.295, width = 0.43, height = 0.32))
  draw_composition()
  popViewport()
  grid.text(c("A", "B", "C"), x = c(0.025, 0.535, 0.025), y = c(0.98, 0.98, 0.48),
            just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
}
png(file.path(out_dir, "Figure2_v29.png"), width = 10.2, height = 13.9358,
    units = "in", res = 300, type = "cairo", bg = "white")
draw_composite()
invisible(dev.off())
cairo_pdf(file.path(out_dir, "Figure2_v29.pdf"), width = 10.2, height = 13.9358)
draw_composite()
invisible(dev.off())
message("GraphPad-style audited Figure 2: 143 groups, 55 barcodes, 363 pairs, 1128 HSPs; no risk classes.")