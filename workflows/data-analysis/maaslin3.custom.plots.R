
# Script: study.R
# Purpose: Generate MaAsLin3 plots for specific variables
# -------------------------------------------------------

# Load required packages
suppressPackageStartupMessages({
  library(maaslin3)
})

# Path to MaAsLin3 output directory
# Use "." if running in same folder as results
output_dir <- "."

# Variables for forest plot (MUST match "metadata value" combined labels, ie. "Timepoint S1" rather than just "Timepoint")
coef_plot_vars <- c("Percussion_Sensitivity Y")

# Variables for heatmaps (optional)
heatmap_vars <- NULL

# Filename for the summary plot PDF
summary_plot_filename <- "S1_Percussion_Sensitivity_plot.pdf"

# Number of top features to plot
summary_plot_first_n <- 30

# Max q-value for significance
max_significance <- 0.1

# Whether to generate individual association plots
plot_associations <- FALSE

# Backup existing summary plot if present
old_plot <- file.path(output_dir, "figures", "summary_plot.pdf")

if (file.exists(old_plot)) {
  timestamp <- format(Sys.time(), "%Y%m%d%H%M")
  backup_plot <- paste0("summary_plot_backup_", timestamp, ".pdf")
  file.copy(old_plot, backup_plot, overwrite = TRUE)
  message("Backup created for existing summary_plot.pdf: ", backup_plot)
}

# Run MaAsLin3's official plotting function
# Run "maaslin_plot_results_from_output" as this utilizes the "study.tsv" output
plots_out <- maaslin_plot_results_from_output(
  output = output_dir,
  metadata = NULL,                  # Not required unless plotting associations
  normalization = "NONE",           # These values are placeholders when re-plotting
  transform = "NONE",
  coef_plot_vars = coef_plot_vars,
  min_beta = 0.2,
  heatmap_vars = heatmap_vars,
  summary_plot_first_n = summary_plot_first_n,
  max_significance = max_significance,
  plot_summary_plot = TRUE,
  plot_associations = plot_associations,
  max_pngs = 30
)

# Rename summary plot to your custom filename
# MaAsLin3 writes summary plot as: "./figures/summary_plot.pdf"

new_plot <- file.path(output_dir, "figures", summary_plot_filename)

if (file.exists(old_plot)) {
  file.rename(old_plot, new_plot)
  message("Custom summary plot saved to: ", new_plot)
} else {
  message("Summary plot was not created (no significant results?).")
  }