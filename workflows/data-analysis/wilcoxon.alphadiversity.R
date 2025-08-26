#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(readr)
})

# CLI argument parser
parser <- ArgumentParser(description = "Alpha diversity analysis: means, Wilcoxon test, FDR correction, boxplots")

parser$add_argument("-i", "--input", required=TRUE, help="Input CSV/TSV file with columns: sample, large_lesion (Y/N), Sobs, Shannon, Chao1")
parser$add_argument("-o", "--output_prefix", required=TRUE, help="Prefix for output files (PDFs and stats table)")
parser$add_argument("--sep", default=",", help="Field separator in input file (default=','). Use '\\t' for tab-delimited.")

args <- parser$parse_args()

# Read input
df <- read_delim(args$input, delim = args$sep, col_types = cols())

# Check required columns
required_cols <- c("sample", "large_lesion", "Sobs", "Shannon", "Chao1")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Function to plot and compute stats
plot_metric <- function(metric, df, output_prefix) {
  formula <- as.formula(paste(metric, "~ large_lesion"))
  test_result <- wilcox.test(formula, data = df)

  p <- ggplot(df, aes(x = large_lesion, y = .data[[metric]], fill = large_lesion)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
    labs(title = paste(metric, "by large_lesion"),
         x = "Large Lesion",
         y = metric) +
    stat_compare_means(method = "wilcox.test", label.y.npc = "top", label = "p.format") +
    theme_minimal(base_size = 14)

  ggsave(filename = paste0(output_prefix, "_", metric, "_boxplot.pdf"),
         plot = p, width = 6, height = 5)

  return(data.frame(
    Metric = metric,
    Group_N = mean(df[[metric]][df$large_lesion == "N"], na.rm = TRUE),
    Group_Y = mean(df[[metric]][df$large_lesion == "Y"], na.rm = TRUE),
    Raw_P = test_result$p.value
  ))
}

# Run for all 3 metrics
results_list <- lapply(c("Sobs", "Shannon", "Chao1"), plot_metric, df = df, output_prefix = args$output_prefix)

results_df <- bind_rows(results_list)

# Apply FDR correction
results_df$Adjusted_P <- p.adjust(results_df$Raw_P, method = "fdr")

# Save stats table
write_csv(results_df, paste0(args$output_prefix, "study.csv"))

cat("✅ Analysis complete!\n")
cat("📁 Outputs:\n")
cat("- PDF plots: ", paste0(args$output_prefix, "_[Sobs|Shannon|Chao1]_boxplot.pdf"), "\n")
cat("- Stats table: ", paste0(args$output_prefix, "study.csv"), "\n")
