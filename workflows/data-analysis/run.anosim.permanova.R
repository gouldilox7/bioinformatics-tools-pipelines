# UNIVERSAL ANOSIM & PERMANOVA ANALYSIS
# Supports .shared, .dist, .csv, and .tsv files

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(tidyverse)
})

# Helper: clean sample names consistently
clean_sample_names <- function(x) {
  x <- gsub("[^[:alnum:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  return(x)
}

# Helper: Convert Mothur long-format .dist to distance matrix
convert_mothur_dist_to_matrix <- function(df) {
  colnames(df) <- c("Sample1", "Sample2", "Distance")
  df$Sample1 <- as.character(df$Sample1)
  df$Sample2 <- as.character(df$Sample2)
  df$Distance <- as.numeric(df$Distance)

  all_samples <- unique(c(df$Sample1, df$Sample2))
  mat <- matrix(0, nrow=length(all_samples), ncol=length(all_samples),
                dimnames=list(all_samples, all_samples))
  for (i in seq_len(nrow(df))) {
    s1 <- df$Sample1[i]
    s2 <- df$Sample2[i]
    d  <- df$Distance[i]
    mat[s1, s2] <- d
    mat[s2, s1] <- d
  }
  as.dist(mat)
}

analyze_dataset <- function(infile, meta_file, output_prefix, nperm = 999) {
  # Read metadata
  meta <- read.delim(meta_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if (!"SampleID" %in% colnames(meta)) stop("Metadata must contain 'SampleID' column.")
  rownames(meta) <- clean_sample_names(meta$SampleID)

  # Load data
  df <- tryCatch({
    if (grepl("\\.csv$", infile, ignore.case=TRUE)) {
      read.csv(infile, check.names=FALSE, stringsAsFactors=FALSE)
    } else {
      read.delim(infile, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
    }
  }, error = function(e) {
    stop("Failed to read file '", infile, "': ", e$message)
  })

  dist_mat <- NULL
  samp_names <- NULL

  if (ncol(df) == 3 && all(sapply(df, is.character) | sapply(df, is.numeric))) {
    message("Detected Mothur .dist file: converting to distance matrix...")
    dist_mat <- convert_mothur_dist_to_matrix(df)
    samp_names <- attr(dist_mat, "Labels")

  } else if (all(c("Group", "numOtus") %in% colnames(df))) {
    message("Detected Mothur .shared file: extracting abundance matrix...")
    rownames(df) <- df$Group
    df <- df[, !(colnames(df) %in% c("label", "Group", "numOtus")), drop=FALSE]
    abund <- df
    samp_names <- clean_sample_names(rownames(abund))
    rownames(abund) <- samp_names
    dist_mat <- vegdist(abund, method="bray")

  } else {
    rownames(df) <- df[[1]]
    df <- df[ , -1, drop=FALSE]
    num_cols <- sapply(df, is.numeric)

    if (all(num_cols) && ncol(df) == nrow(df)) {
      message("Detected square numeric matrix: treating as precomputed distance matrix...")
      dist_mat <- as.dist(as.matrix(df))
      samp_names <- rownames(df)
    } else {
      message("Detected abundance matrix: computing Bray-Curtis distance...")
      abund <- df[, num_cols, drop=FALSE]
      if (ncol(abund) < 1) stop("No numeric abundance columns found in ", infile)
      abund_t <- t(abund)
      samp_names <- clean_sample_names(rownames(abund_t))
      rownames(abund_t) <- samp_names
      dist_mat <- vegdist(abund_t, method="bray")
    }
  }

  # Apply cleaned names to distance matrix
  if (inherits(dist_mat, "dist")) {
    attr(dist_mat, "Labels") <- samp_names
  }

  # Subset metadata and distance matrix
  rownames(meta) <- clean_sample_names(rownames(meta))
  common <- intersect(attr(dist_mat, "Labels"), rownames(meta))

  if (length(common) < 2) {
    warning("Skipping file '", infile, "': fewer than 2 samples in common with metadata.")
    return(NULL)
  }

  dist_mat <- as.dist(as.matrix(dist_mat)[common, common])
  meta <- meta[common, , drop=FALSE]

  # Run ANOSIM and PERMANOVA
  vars <- setdiff(colnames(meta), "SampleID")
  results <- list()

  for (var in vars) {
    grp <- meta[[var]]
    if (all(is.na(grp)) || length(unique(grp[!is.na(grp)])) < 2) {
      warning("Skipping variable '", var, "': fewer than 2 groups or all NA.")
      next
    }
    meta$.group <- factor(grp)

    a <- anosim(dist_mat, meta$.group, permutations=nperm)
    p <- adonis2(dist_mat ~ .group, data=meta, permutations=nperm, method="bray")[1, ]

    results[[var]] <- tibble(
      Variable     = var,
      Test         = c("ANOSIM", "PERMANOVA"),
      R_statistic  = c(a$statistic, NA_real_),
      P_value      = c(a$signif, p$`Pr(>F)`),
      Df           = c(NA_integer_, p$Df),
      SumOfSqs     = c(NA_real_, p$SumOfSqs),
      R2           = c(NA_real_, p$R2),
      F            = c(NA_real_, p$F),
      Permutations = c(nperm, nperm)
    )
  }

  out_df <- bind_rows(results)
  out_file <- paste0(output_prefix, "study.tsv")
  write.table(out_df, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  message("Wrote results: ", out_file)
}

# CLI Options
option_list <- list(
  make_option(c("-i","--input"), type="character", help="Comma-separated input files", metavar="files"),
  make_option(c("-m","--meta"),  type="character", help="Metadata file (TSV with SampleID)", metavar="file"),
  make_option(c("-o","--out"),   type="character", default="results", help="Output prefix [default %default]", metavar="prefix"),
  make_option(c("-n","--nperm"), type="integer", default=999, help="Permutations [default %default]", metavar="num")
)

parser <- OptionParser(option_list=option_list)
opts <- parse_args(parser)
if (is.null(opts$input) || is.null(opts$meta)) {
  print_help(parser)
  stop("--input and --meta are required.", call.=FALSE)
}

files <- strsplit(opts$input, ",")[[1]]
for (f in files) {
  message("Analyzing file: ", f)
  prefix <- paste(opts$out, tools::file_path_sans_ext(basename(f)), sep="_")
  analyze_dataset(f, opts$meta, prefix, opts$nperm)
}
