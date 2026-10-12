#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(bitmapType = "cairo")

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[[1]] else getwd()
indir <- file.path(base, "19_local_genealogy", "results")
outdir <- file.path(base, "19_local_genealogy", "plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

metrics_file <- file.path(indir, "all650_diploid_distance_metrics.tsv")
tests_file <- file.path(indir, "genealogy_class_permutation_tests.tsv")
stopifnot(file.exists(metrics_file), file.exists(tests_file))

dat <- fread(metrics_file, na.strings = c("NA", "NaN", ""))
tests <- fread(tests_file, na.strings = c("NA", "NaN", ""))

model_levels <- c(
  "balancing_selection",
  "recurrent_selection",
  "allopatric_selection",
  "divergence_with_gene_flow",
  "unclassified"
)
model_labels <- c(
  divergence_with_gene_flow = "Div. + GF",
  allopatric_selection = "Allopatric",
  recurrent_selection = "Recurrent",
  balancing_selection = "Balancing",
  unclassified = "Unclass."
)
model_palette <- c(
  balancing_selection = "#7B3294",
  recurrent_selection = "#0072B2",
  allopatric_selection = "#E69F00",
  divergence_with_gene_flow = "#D73027",
  unclassified = "grey78"
)

long <- melt(
  dat,
  id.vars = c("locus_id", "model_class"),
  measure.vars = c("relative_mean_depth_median", "nearest_neighbor_mixing_median"),
  variable.name = "metric_code",
  value.name = "value",
  variable.factor = FALSE
)
long <- long[is.finite(value) & model_class %chin% model_levels]
long[, model_class := factor(model_class, levels = model_levels)]
long[, metric := factor(
  metric_code,
  levels = c("relative_mean_depth_median", "nearest_neighbor_mixing_median"),
  labels = c(
    "Reference-normalized mean pairwise depth",
    "PAP-PSA nearest-neighbor mixing"
  )
)]

# Counts are identical between panels for the current complete-case data, but
# labels are calculated from the input rather than hard-coded.
counts <- unique(long[, .(locus_id, model_class)])[, .N, by = model_class]
count_lookup <- setNames(counts$N, as.character(counts$model_class))
x_labels <- setNames(
  vapply(model_levels, function(x) {
    sprintf("%s\n(n = %d)", model_labels[[x]], count_lookup[[x]])
  }, character(1)),
  model_levels
)

# A one-locus class cannot support a density estimate. Retain it in the
# boxplot layer and omit only its undefined violin density.
violin_dat <- long[, if (.N >= 2L) .SD, by = .(metric, model_class)]

get_q <- function(metric_name) {
  z <- tests[
    analysis == "all_650_phase_free" &
      metric == metric_name &
      comparison == "balancing_vs_unclassified"
  ]
  if (nrow(z) != 1L) return(NA_real_)
  z$q_BH_within_analysis[[1]]
}

depth_q <- get_q("relative_mean_depth_median")
mixing_q <- get_q("nearest_neighbor_mixing_median")
subtitle <- sprintf(
  "Balancing vs unclassified permutation tests: depth BH Q = %.4f; mixing BH Q = %.4f",
  depth_q, mixing_q
)

p <- ggplot(long, aes(model_class, value, fill = model_class, color = model_class)) +
  geom_violin(
    data = violin_dat,
    scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.58,
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.75,
    na.rm = TRUE
  ) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = model_palette, drop = FALSE, guide = "none") +
  scale_color_manual(values = model_palette, drop = FALSE, guide = "none") +
  scale_x_discrete(labels = x_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = "Locus statistic",
    title = "PAP-PSA local-genealogy distributions across 650 RAD loci",
    subtitle = subtitle
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11.5),
    axis.title.y = element_text(size = 13),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10.5, angle = 0, hjust = 0.5, lineheight = 0.9),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    plot.margin = margin(8, 10, 8, 8)
  )

preferred_pdf <- file.path(outdir, "pap_psa_all650_genealogy_comparison_main.pdf")
preferred_png <- file.path(outdir, "pap_psa_all650_genealogy_comparison_main.png")
r_pdf <- file.path(outdir, "pap_psa_all650_genealogy_comparison_main_R.pdf")
r_png <- file.path(outdir, "pap_psa_all650_genealogy_comparison_main_R.png")

ggsave(preferred_pdf, p, width = 8.5, height = 4.8, useDingbats = FALSE, bg = "white")
ggsave(preferred_png, p, width = 8.5, height = 4.8, dpi = 240, bg = "white")
ggsave(r_pdf, p, width = 8.5, height = 4.8, useDingbats = FALSE, bg = "white")
ggsave(r_png, p, width = 8.5, height = 4.8, dpi = 240, bg = "white")

summary_table <- long[, .(
  n = .N,
  median = median(value),
  mean = mean(value),
  sd = sd(value),
  q025 = quantile(value, 0.025),
  q975 = quantile(value, 0.975)
), by = .(metric, model_class)]
fwrite(
  summary_table,
  file.path(indir, "all650_genealogy_main_figure_summary.tsv"),
  sep = "\t"
)
