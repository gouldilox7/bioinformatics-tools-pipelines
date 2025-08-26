#!/usr/bin/env Rscript
# ======================================================================
# run.anosim.permanova.R  —  One‑stop CLI for beta‑diversity significance testing
# ======================================================================
# This script inspects your input distance/abundance table + metadata,
# runs ANOSIM, PERMANOVA, and PERMDISP (overall + pairwise, with optional nesting),
# and writes clean TSVs:
#   <prefix>_ANOSIM_PERMANOVA_PERMDISP.tsv
#   <prefix>_ANOSIM_PERMANOVA_PERMDISP_pairwise.tsv
#
# SUPPORTED INPUTS
# ----------------
# 1) Mothur .dist (long format; 3 columns: Sample1, Sample2, Distance)
# 2) Mothur .shared (rows = samples, columns = OTUs; Bray–Curtis computed)
# 3) Square numeric distance matrix (row/col names are samples)
# 4) Abundance table (first column feature ID; remaining columns sample IDs; Bray–Curtis computed)
#
# TESTS
# -----
# • ANOSIM (vegan::anosim): rank‑based between/within distances → R statistic + p
# • PERMANOVA (vegan::adonis2): variance partitioning on distances → F, R2, p
# • PERMDISP (vegan::betadisper + permutest): dispersion (homogeneity of multivariate dispersion)
#
# DESIGN
# ------
# • GLOBAL (overall) tests for each non‑numeric metadata column (excluding SampleID)
# • PAIRWISE tests across all level pairs for each factor (guarded by min-per-group)
# • NESTED tests (optional): run a TEST variable within each level of BY variable
#   Using CLI like:  --nested "Group|Generation,Generation|Group"
#
# OUTPUT COLUMNS (overall)
# ------------------------
# Variable, ByVar, ByLevel, Contrast, Test, R_statistic, P_value, Df, SumOfSqs, R2, F, Permutations
#
# OUTPUT COLUMNS (pairwise)
# -------------------------
# Variable, ByVar, ByLevel, Contrast, Group1, Group2, N1, N2, Test, R_statistic, P_value, Df, SumOfSqs, R2, F, Permutations
#
# CONVENTIONS
# -----------
# • Metadata must contain a 'SampleID' column (case‑insensitive is OK; will be normalized)
# • Sample names are cleaned to safe identifiers: non‑alnum → '_', collapse runs, trim edges
# • Distance computation uses Bray–Curtis for abundance‑type inputs
# • All outputs are deterministically ordered (by p‑value ascending within Test)
#
# EXAMPLES
# --------
#  Basic (one input file):
#    Rscript run.anosim.permanova.R --input table.tsv --meta meta.tsv --out results
#
#  Multiple inputs (comma‑separated):
#    Rscript run.anosim.permanova.R -i "A.tsv,B.tsv" -m meta.tsv -o results
#
#  With nesting (test Group within each Generation, and vice‑versa):
#    Rscript run.anosim.permanova.R -i table.tsv -m meta.tsv --nested "Group|Generation,Generation|Group"
#
#  Custom permutations and diagnostic log:
#    Rscript run.anosim.permanova.R -i table.tsv -m meta.tsv -n 4999 --log beta_tests.log --diag
#
# ======================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(tidyverse)
})

# ---------- Dependency helper (consistent with alpha script style) -------
ensure_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s. Install with install.packages('%s')",
                 pkg, why, pkg), call. = FALSE)
  }
}
ensure_pkg("optparse", "CLI parsing")
ensure_pkg("vegan", "ANOSIM/PERMANOVA/PERMDISP")
ensure_pkg("tibble", "tidy tibble outputs")
ensure_pkg("dplyr", "data manipulation")
ensure_pkg("tidyr", "data reshaping")
ensure_pkg("readr", "robust I/O")

# ---------- Utility & logging -------------------------------------------
basename_noext <- function(p) tools::file_path_sans_ext(basename(p))

pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error=function(e) NA_character_)

# Set by CLI if --log provided (see below)
.log_path <- NULL
log_line <- function(...) {
  if (!is.null(.log_path) && nzchar(.log_path)) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste(ts, "-", paste(..., collapse = " ")), "\n", file = .log_path, append = TRUE)
  }
}

# ---------- Sample name cleaner (matches alpha script ethos) -------------
clean_sample_names <- function(x) {
  x <- gsub("[^[:alnum:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---------- Mothur .dist → 'dist' object --------------------------------
convert_mothur_dist_to_matrix <- function(df) {
  colnames(df) <- c("Sample1", "Sample2", "Distance")
  df$Sample1 <- clean_sample_names(as.character(df$Sample1))
  df$Sample2 <- clean_sample_names(as.character(df$Sample2))
  df$Distance <- as.numeric(df$Distance)

  all_samples <- sort(unique(c(df$Sample1, df$Sample2)))
  mat <- matrix(0, nrow=length(all_samples), ncol=length(all_samples),
                dimnames=list(all_samples, all_samples))
  for (i in seq_len(nrow(df))) {
    s1 <- df$Sample1[i]; s2 <- df$Sample2[i]; d <- df$Distance[i]
    mat[s1, s2] <- d; mat[s2, s1] <- d
  }
  as.dist(mat)
}

# ---------- Pairwise helper (carries nested context columns) -------------
pairwise_beta <- function(dist_full, meta_df, var, nperm=999, min_per_group=2,
                          by_var=NA_character_, by_level=NA_character_, contrast=NA_character_) {
  # Prepare grouping factor
  grp <- meta_df[[var]]
  names(grp) <- rownames(meta_df)
  keep <- !is.na(grp)
  grp <- droplevels(factor(grp[keep]))
  if (nlevels(grp) < 2) return(tibble())

  # Align to distance labels
  labs <- attr(dist_full, "Labels")
  common <- intersect(names(grp), labs)
  if (length(common) < 2) return(tibble())
  dist_full <- as.dist(as.matrix(dist_full)[common, common])
  grp <- grp[common]

  # All level pairs
  lvl <- levels(grp)
  combs <- combn(lvl, 2, simplify=FALSE)

  res_list <- list()
  for (ab in combs) {
    g1 <- ab[1]; g2 <- ab[2]
    idx <- grp %in% c(g1, g2)
    sub_grp <- droplevels(grp[idx])
    if (nlevels(sub_grp) < 2) next
    tab <- table(sub_grp)
    if (any(tab < min_per_group)) next

    sub_labels <- names(sub_grp)
    sub_dist <- as.dist(as.matrix(dist_full)[sub_labels, sub_labels])
    sub_meta <- tibble(.group = sub_grp, rowname = sub_labels) |>
      column_to_rownames("rowname")

    # ANOSIM
    a <- try(anosim(sub_dist, sub_meta$.group, permutations=nperm), silent=TRUE)
    # PERMANOVA (method ignored when 'dist' provided; kept for clarity)
    p <- try(adonis2(sub_dist ~ .group, data=sub_meta, permutations=nperm, method="bray"),
             silent=TRUE)
    # PERMDISP
    permdisp_F <- NA_real_; permdisp_p <- NA_real_; permdisp_df <- NA_integer_
    bd <- try(betadisper(sub_dist, sub_meta$.group), silent=TRUE)
    if (!inherits(bd, "try-error")) {
      pt <- try(permutest(bd, permutations=nperm), silent=TRUE)
      if (!inherits(pt, "try-error")) {
        tab_row <- tryCatch(as.data.frame(pt$tab)[1, , drop=FALSE], error=function(e) NULL)
        if (!is.null(tab_row)) {
          permdisp_F  <- suppressWarnings(as.numeric(tab_row$F))
          pcol <- intersect(colnames(tab_row), c("Pr(>F)","Pr(>F)1","Pr(>F)2","Pr(>F)3"))
          if (length(pcol) == 0) pcol <- colnames(tab_row)[grepl("Pr\\(>F\\)", colnames(tab_row))]
          if (length(pcol)) permdisp_p <- suppressWarnings(as.numeric(tab_row[[pcol[1]]]))
          permdisp_df <- suppressWarnings(as.integer(tab_row$Df))
        }
      }
    }

    one_pair <- tibble(
      Variable     = var,
      ByVar        = by_var,
      ByLevel      = by_level,
      Contrast     = contrast,
      Group1       = g1,
      Group2       = g2,
      N1           = unname(tab[g1]),
      N2           = unname(tab[g2]),
      Test         = c("ANOSIM", "PERMANOVA", "PERMDISP"),
      R_statistic  = c(if (inherits(a,"anosim")) a$statistic else NA_real_, NA_real_, NA_real_),
      P_value      = c(if (inherits(a,"anosim")) a$signif else NA_real_,
                       if (inherits(p,"anova")) p$`Pr(>F)`[1] else NA_real_,
                       permdisp_p),
      Df           = c(NA_integer_,
                       if (inherits(p,"anova")) p$Df[1] else NA_integer_,
                       permdisp_df),
      SumOfSqs     = c(NA_real_,
                       if (inherits(p,"anova")) p$SumOfSqs[1] else NA_real_,
                       NA_real_),
      R2           = c(NA_real_,
                       if (inherits(p,"anova")) p$R2[1] else NA_real_,
                       NA_real_),
      F            = c(NA_real_,
                       if (inherits(p,"anova")) p$F[1] else NA_real_,
                       permdisp_F),
      Permutations = c(nperm, nperm, nperm)
    )
    res_list[[length(res_list)+1]] <- one_pair
  }
  bind_rows(res_list)
}

# ---------- Core block runner (overall + delegates pairwise) -------------
run_overall_block <- function(dist_mat, meta, var, nperm=999,
                              by_var=NA_character_, by_level=NA_character_, contrast=NA_character_) {
  grp <- meta[[var]]
  if (all(is.na(grp)) || length(unique(grp[!is.na(grp)])) < 2) return(list(overall=NULL, pw=NULL))
  meta$.group <- factor(grp)

  # ANOSIM
  a <- try(anosim(dist_mat, meta$.group, permutations=nperm), silent=TRUE)
  # PERMANOVA
  p <- try(adonis2(dist_mat ~ .group, data=meta, permutations=nperm, method="bray"), silent=TRUE)
  # PERMDISP
  permdisp_F <- NA_real_; permdisp_p <- NA_real_; permdisp_df <- NA_integer_
  bd <- try(betadisper(dist_mat, meta$.group), silent=TRUE)
  if (!inherits(bd, "try-error")) {
    pt <- try(permutest(bd, permutations=nperm), silent=TRUE)
    if (!inherits(pt, "try-error")) {
      tab_row <- tryCatch(as.data.frame(pt$tab)[1, , drop=FALSE], error=function(e) NULL)
      if (!is.null(tab_row)) {
        permdisp_F  <- suppressWarnings(as.numeric(tab_row$F))
        pcol <- intersect(colnames(tab_row), c("Pr(>F)","Pr(>F)1","Pr(>F)2","Pr(>F)3"))
        if (length(pcol) == 0) pcol <- colnames(tab_row)[grepl("Pr\\(>F\\)", colnames(tab_row))]
        if (length(pcol)) permdisp_p <- suppressWarnings(as.numeric(tab_row[[pcol[1]]]))
        permdisp_df <- suppressWarnings(as.integer(tab_row$Df))
      }
    }
  }

  overall <- tibble(
    Variable     = var,
    ByVar        = by_var,
    ByLevel      = by_level,
    Contrast     = contrast,
    Test         = c("ANOSIM", "PERMANOVA", "PERMDISP"),
    R_statistic  = c(if (inherits(a,"anosim")) a$statistic else NA_real_, NA_real_, NA_real_),
    P_value      = c(if (inherits(a,"anosim")) a$signif else NA_real_,
                     if (inherits(p,"anova")) p$`Pr(>F)`[1] else NA_real_,
                     permdisp_p),
    Df           = c(NA_integer_,
                     if (inherits(p,"anova")) p$Df[1] else NA_integer_,
                     permdisp_df),
    SumOfSqs     = c(NA_real_,
                     if (inherits(p,"anova")) p$SumOfSqs[1] else NA_real_,
                     NA_real_),
    R2           = c(NA_real_,
                     if (inherits(p,"anova")) p$R2[1] else NA_real_,
                     NA_real_),
    F            = c(NA_real_,
                     if (inherits(p,"anova")) p$F[1] else NA_real_,
                     permdisp_F),
    Permutations = c(nperm, nperm, nperm)
  )

  pw <- pairwise_beta(dist_mat, meta, var, nperm=nperm,
                      by_var=by_var, by_level=by_level, contrast=contrast)

  list(overall=overall, pw=pw)
}

# ---------- Dataset analyzer (per input file) ----------------------------
analyze_dataset <- function(infile, meta_file, output_prefix, nperm = 999, nested_specs=NULL,
                            min_per_group = 2, diag = FALSE) {

  log_line("ANALYZE", infile, "meta=", meta_file, "nperm=", nperm, "diag=", diag)

  # --- Read metadata (normalize SampleID name) ---
  meta <- readr::read_tsv(meta_file, show_col_types = FALSE)
  nm_low <- tolower(names(meta))
  if (!"sampleid" %in% nm_low) stop("Metadata must contain a 'SampleID' column (case-insensitive).")
  names(meta)[nm_low=="sampleid"] <- "SampleID"
  meta$SampleID <- clean_sample_names(meta$SampleID)
  rownames(meta) <- meta$SampleID

  # --- Load primary table ---
  df <- tryCatch({
    ext <- tolower(tools::file_ext(infile))
    if (ext %in% c("csv")) {
      readr::read_csv(infile, show_col_types = FALSE)
    } else {
      readr::read_tsv(infile, show_col_types = FALSE)
    }
  }, error = function(e) stop("Failed to read file '", infile, "': ", e$message))

  dist_mat <- NULL
  samp_names <- NULL

  # --- Detect input type & compute/convert distances ---
  if (ncol(df) == 3 && all(sapply(df, function(x) is.character(x) || is.numeric(x)))) {
    message("Detected Mothur .dist file: converting to distance matrix...")
    log_line("Input detected: mothur .dist")
    dist_mat <- convert_mothur_dist_to_matrix(df)
    samp_names <- attr(dist_mat, "Labels")

  } else if (all(c("Group", "numOtus") %in% colnames(df))) {
    message("Detected Mothur .shared file: extracting abundance matrix...")
    log_line("Input detected: mothur .shared")
    rownames(df) <- df$Group
    df <- df[, !(colnames(df) %in% c("label", "Group", "numOtus")), drop=FALSE]
    abund <- df
    samp_names <- clean_sample_names(rownames(abund))
    rownames(abund) <- samp_names
    dist_mat <- vegdist(abund, method="bray")

  } else {
    # If first column looks like sample IDs (square matrix) or feature IDs (abundance)
    rownames(df) <- clean_sample_names(df[[1]])
    df <- df[, -1, drop=FALSE]
    num_cols <- sapply(df, is.numeric)

    if (all(num_cols) && ncol(df) == nrow(df)) {
      message("Detected square numeric matrix: treating as precomputed distance matrix...")
      log_line("Input detected: square distance matrix")
      dist_mat <- as.dist(as.matrix(df))
      samp_names <- rownames(df)
    } else {
      message("Detected abundance matrix: computing Bray-Curtis distance...")
      log_line("Input detected: abundance table (Bray–Curtis)")
      abund <- df[, num_cols, drop=FALSE]
      if (ncol(abund) < 1) stop("No numeric abundance columns found in ", infile)
      abund_t <- t(abund)
      samp_names <- clean_sample_names(rownames(abund_t))
      rownames(abund_t) <- samp_names
      dist_mat <- vegdist(abund_t, method="bray")
    }
  }

  if (inherits(dist_mat, "dist")) attr(dist_mat, "Labels") <- clean_sample_names(samp_names)

  # --- Align with metadata ---
  common <- intersect(attr(dist_mat, "Labels"), rownames(meta))
  if (length(common) < 2) {
    warning("Skipping file '", infile, "': fewer than 2 samples in common with metadata.")
    log_line("WARN: fewer than 2 samples overlap; skipping", infile)
    return(NULL)
  }
  dist_mat <- as.dist(as.matrix(dist_mat)[common, common])
  meta <- meta[common, , drop=FALSE]

  # --- Variables to test (non‑numeric columns except SampleID) ---
  vars <- setdiff(colnames(meta), "SampleID")
  # Proactive drop numeric columns to avoid accidental testing of covariates
  vars <- vars[!sapply(meta[vars], is.numeric)]

  results_overall <- list()
  results_pairwise <- list()

  # --- GLOBAL analyses ---
  for (var in vars) {
    block <- run_overall_block(dist_mat, meta, var, nperm=nperm,
                               by_var=NA_character_, by_level=NA_character_, contrast=NA_character_)
    if (!is.null(block$overall)) results_overall[[paste0("global_",var)]] <- block$overall
    if (!is.null(block$pw) && nrow(block$pw)) results_pairwise[[paste0("global_",var)]] <- block$pw
  }

  # --- NESTED analyses (if requested) ---
  if (!is.null(nested_specs) && length(nested_specs)) {
    for (spec in nested_specs) {
      parts <- strsplit(spec, "\\|")[[1]]
      if (length(parts) != 2) { warning("Ignoring malformed --nested spec: ", spec); next }
      test_var <- parts[1]; by_var <- parts[2]

      if (!test_var %in% colnames(meta)) { warning("Test var '", test_var, "' not in metadata; skipping."); next }
      if (!by_var %in% colnames(meta))   { warning("By var '", by_var,   "' not in metadata; skipping."); next }
      lvls <- sort(unique(meta[[by_var]]))

      for (lv in lvls) {
        sub_idx <- which(meta[[by_var]] == lv)
        if (length(sub_idx) < 2) next
        meta_sub <- meta[sub_idx, , drop=FALSE]

        # Align distance subset
        labs <- attr(dist_mat, "Labels")
        common_sub <- intersect(rownames(meta_sub), labs)
        if (length(common_sub) < 2) next
        dist_sub <- as.dist(as.matrix(dist_mat)[common_sub, common_sub])
        meta_sub <- meta_sub[common_sub, , drop=FALSE]

        # Run block for test_var within this level of by_var
        block <- run_overall_block(dist_sub, meta_sub, test_var, nperm=nperm,
                                   by_var=by_var, by_level=as.character(lv), contrast=test_var)

        key <- paste0("nested_", test_var, "_within_", by_var, "_", lv)
        if (!is.null(block$overall)) results_overall[[key]] <- block$overall
        if (!is.null(block$pw) && nrow(block$pw)) results_pairwise[[key]] <- block$pw
      }
    }
  }

  # --- Bind & sort (deterministic, p‑value ascending within Test) ---
  out_overall <- bind_rows(results_overall) %>%
    arrange(Test, is.na(P_value), P_value)
  out_pw <- bind_rows(results_pairwise) %>%
    arrange(Test, is.na(P_value), P_value)

  # --- Optional diagnostics: counts per factor level (global only) ---
  if (isTRUE(diag)) {
    diag_rows <- list()
    for (var in vars) {
      tab <- table(meta[[var]])
      if (length(tab)) {
        diag_rows[[length(diag_rows)+1]] <- tibble(Variable = var,
                                                   Level = names(tab),
                                                   N = as.integer(tab))
      }
    }
    if (length(diag_rows)) {
      readr::write_tsv(bind_rows(diag_rows), paste0(output_prefix, "_diag_counts.tsv"))
      log_line("WROTE", paste0(output_prefix, "_diag_counts.tsv"))
    }
  }

  # --- Write outputs ---
  out_file_overall <- paste0(output_prefix, "_ANOSIM_PERMANOVA_PERMDISP.tsv")
  readr::write_tsv(out_overall, out_file_overall, na = "")
  message("Wrote results: ", out_file_overall)
  log_line("WROTE", out_file_overall, "rows=", nrow(out_overall))

  out_file_pw <- paste0(output_prefix, "_ANOSIM_PERMANOVA_PERMDISP_pairwise.tsv")
  if (nrow(out_pw)) {
    readr::write_tsv(out_pw, out_file_pw, na = "")
    message("Wrote pairwise results: ", out_file_pw)
    log_line("WROTE", out_file_pw, "rows=", nrow(out_pw))
  } else {
    message("No pairwise results produced (insufficient levels or sample sizes).")
    log_line("INFO: no pairwise results")
  }
}

# ---------- CLI Options (documented like alpha script) -------------------
option_list <- list(
  make_option(c("-i","--input"), type="character",
              help="Comma-separated input files (required). Examples: 'A.tsv,B.tsv' or 'table.tsv'"),
  make_option(c("-m","--meta"),  type="character",
              help="Metadata file (TSV/CSV) with 'SampleID' column (required)."),
  make_option(c("-o","--out"),   type="character", default="results",
              help="Output prefix [default: %default]"),
  make_option(c("-n","--nperm"), type="integer",   default=999,
              help="Number of permutations [default: %default]"),
  make_option(c("--nested"),     type="character", default=NULL,
              help="Directional nested specs as TEST|BY (comma-separated). Example: 'Group|Generation,Generation|Group'"),
  make_option(c("--min_per_group"), type="integer", default=2,
              help="Minimum samples per level to run pairwise tests [default: %default]"),
  make_option(c("--diag"), action="store_true", default=FALSE,
              help="Write simple diagnostic TSVs (e.g., counts per global factor)."),
  make_option(c("--log"), type="character", default=NULL,
              help="Path to diagnostic log file. If supplied, progress & environment are logged.")
)

parser <- OptionParser(option_list=option_list)
opts <- parse_args(parser)

# Required checks (clear errors like alpha script)
if (is.null(opts$input) || is.null(opts$meta)) {
  print_help(parser)
  stop("--input and --meta are required.", call.=FALSE)
}

# Initialize logging
start_time <- Sys.time()
if (!is.null(opts$log) && nzchar(opts$log)) {
  .log_path <<- opts$log
  log_line("START run.anosim.permanova.R")
  log_line("R:", R.version.string)
  log_line("Packages:", paste0(
    "vegan=", pkg_ver("vegan"), ", dplyr=", pkg_ver("dplyr"),
    ", tidyr=", pkg_ver("tidyr"), ", readr=", pkg_ver("readr"),
    ", optparse=", pkg_ver("optparse")))
  log_line("Options parsed.")
}

# Parse nested specs string → vector
nested_specs <- NULL
if (!is.null(opts$nested) && nzchar(opts$nested)) {
  nested_specs <- strsplit(opts$nested, "\\s*,\\s*")[[1]]
  nested_specs <- nested_specs[nzchar(nested_specs)]
}

# Fan out over inputs
files <- strsplit(opts$input, ",")[[1]]
for (f in files) {
  f <- trimws(f)
  if (!nzchar(f)) next
  message("Analyzing file: ", f)
  log_line("FILE", f)
  prefix <- paste(opts$out, basename_noext(f), sep="_")
  analyze_dataset(f, opts$meta, prefix, opts$nperm, nested_specs,
                  min_per_group = opts$min_per_group, diag = isTRUE(opts$diag))
}

if (!is.null(.log_path) && nzchar(.log_path)) {
  dt <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  log_line(sprintf("DONE in %.2f sec", dt))
}
message("✅ Done.")
