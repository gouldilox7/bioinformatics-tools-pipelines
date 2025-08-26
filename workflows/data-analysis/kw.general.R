
# Load libraries
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))

# CLI options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Feature table (tab-delimited) [required]"),
  make_option(c("-m", "--metadata"), type = "character", help = "Metadata file with sample IDs and grouping [required]"),
  make_option(c("-c", "--group_column"), type = "character", help = "Column in metadata that defines the groups [required]"),
  make_option(c("-g", "--group1"), type = "character", help = "First group label [required]"),
  make_option(c("-G", "--group2"), type = "character", help = "Second group label [required]"),
  make_option(c("-o", "--output"), type = "character", default = "study.csv", help = "Output file [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate input
required <- c("input", "metadata", "group_column", "group1", "group2")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", "), "\nUse -h for help."))
}

# Load data
feature_table <- read_tsv(opt$input)
metadata <- read_tsv(opt$metadata)

# Sample names = all columns except first
sample_ids <- colnames(feature_table)[-1]

# Filter metadata
metadata <- metadata %>%
  filter(!!sym(opt$group_column) %in% c(opt$group1, opt$group2)) %>%
  filter(sample %in% sample_ids)

# Get sample IDs by group
group1_samples <- metadata %>%
  filter(!!sym(opt$group_column) == opt$group1) %>%
  pull(sample)

group2_samples <- metadata %>%
  filter(!!sym(opt$group_column) == opt$group2) %>%
  pull(sample)

if (length(group1_samples) == 0 | length(group2_samples) == 0) {
  stop("❌ One or both groups have no samples after filtering.")
}

# Perform Kruskal-Wallis tests
results_list <- vector("list", nrow(feature_table))

for (i in seq_len(nrow(feature_table))) {
  row_data <- feature_table[i, ]
  feature_id <- row_data[[1]]
  group1_vals <- as.numeric(row_data[, group1_samples, drop = FALSE])
  group2_vals <- as.numeric(row_data[, group2_samples, drop = FALSE])

  # Only test if variance exists
  if (length(unique(c(group1_vals, group2_vals))) > 1) {
    test <- tryCatch(kruskal.test(list(group1_vals, group2_vals)), error = function(e) NULL)
    if (!is.null(test)) {
      results_list[[i]] <- tibble(
        Feature = feature_id,
        Group1_Avg = mean(group1_vals, na.rm = TRUE),
        Group2_Avg = mean(group2_vals, na.rm = TRUE),
        Statistic = as.numeric(test$statistic),
        P_value = test$p.value
      )
    }
  }
}

# Combine and finalize
results_df <- bind_rows(results_list)

if (nrow(results_df) > 0) {
  results_df <- results_df %>%
    mutate(FDR_adjusted_p = p.adjust(P_value, method = "BH")) %>%
    arrange(P_value) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Feature, Statistic, P_value, FDR_adjusted_p, Group1_Avg, Group2_Avg)

  write_csv(results_df, opt$output)
  cat("✅ Kruskal-Wallis test complete.\n📄 Results saved to:", opt$output, "\n")
} else {
  cat("⚠️ No valid rows with variance between groups. No output file written.\n")
}
