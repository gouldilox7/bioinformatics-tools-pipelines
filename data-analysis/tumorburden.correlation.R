#!/usr/bin/env Rscript

#===========================
# tumorburden.correlation.R
#===========================
# Correlate tumor burden (MouseID-level, longitudinal) to metagenomic data (sample-level).
#
# Outputs (default --outdir outputs):
#  - outputs/intermediate/master_metadata.tsv
#  - outputs/intermediate/master_alpha.tsv
#  - outputs/intermediate/master_genus_rel.tsv
#  - outputs/endpoint/alpha_vs_tumor.tsv          (Spearman rho/p/q)
#  - outputs/endpoint/alpha_lmm_vs_tumor.tsv      (if Cohort present)
#  - outputs/endpoint/envfit_tumor.tsv            (vector fit r^2, p)
#  - outputs/endpoint/mantel_tumor.tsv            (Mantel r, p, permutations)
#  - outputs/maaslin3_tumor/                      (MaAsLin3 outputs)
#
# Example:
#  Rscript tumorburden.correlation.R \
#    --metadata_files 1st.fmt.metadata.txt \
#    --metadata_files 1st.second.metadata.txt \
#    --metadata_files 1st.third.metadata.txt \
#    --metadata_files 3rd.metadata.txt \
#    --alpha_files 1st.fmt.alphametrics.txt \
#    --alpha_files 1st.second.alphametrics.txt \
#    --alpha_files 1st.third.alphametrics.txt \
#    --alpha_files 3rd.alphametrics.txt \
#    --genus_files 1st.fmt.genustable.txt \
#    --genus_files 1st.second.genustable.txt \
#    --genus_files 1st.third.genustable.txt \
#    --genus_files 3rd.genustable.txt \
#    --outdir outputs
#
suppressPackageStartupMessages({
  library(optparse)     # CRAN: CLI parsing
  library(readr)        # CRAN: TSV/CSV IO
  library(readxl)       # CRAN: xlsx support (optional)
  library(dplyr)        # CRAN: data wrangling
  library(tidyr)        # CRAN: pivoting
  library(stringr)      # CRAN: strings
  library(purrr)        # CRAN: mapping
  library(tibble)       # CRAN: rownames<->column helpers
  library(vegan)        # CRAN: envfit, mantel, vegdist, capscale
  library(lme4)         # CRAN: lmer mixed models
  library(broom.mixed)  # CRAN: tidy lmer results
  library(maaslin3)     # Bioconductor: MaAsLin3 associations
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

die <- function(msg) { cat(sprintf("\\n[FATAL] %s\\n\\n", msg)); quit(status = 1) }
inform <- function(msg) { cat(sprintf("[INFO] %s\\n", msg)) }
warnf <- function(msg) { cat(sprintf("[WARN] %s\\n", msg)) }

# ---------- CLI ----------
opt_list <- list(
  make_option(c("--tumor_file"), type="character", default="tumor.metadata.txt",
              help="Path to tumor burden file (.tsv/.txt/.csv or .xlsx). [default %default]"),
  make_option(c("--metadata_files"), type="character", action="append",
              help="One or more cohort metadata files. Repeat flag to add multiple. Required."),
  make_option(c("--alpha_files"), type="character", action="append", default=NULL,
              help="One or more alpha metrics files. Repeat flag to add multiple. Optional."),
  make_option(c("--genus_files"), type="character", action="append",
              help="One or more genus tables (genera rows, sample columns). Repeat to add multiple. Required."),
  make_option(c("--outdir"), type="character", default="outputs",
              help="Output directory [default %default]"),
  make_option(c("--endpoint_regex"), type="character", default="(?i)(sac|end|terminal)",
              help="Regex to detect endpoint stool sample in Microbiome Timepoint [default %default]"),
  make_option(c("--permutations"), type="integer", default=9999,
              help="Permutations for envfit/mantel [default %default]"),
  make_option(c("--seed"), type="integer", default=42,
              help="Random seed for reproducibility [default %default]"),
  make_option(c("--run_alpha"), type="logical", default=TRUE,
              help="Run alpha ~ tumor correlation [default %default]"),
  make_option(c("--run_beta"), type="logical", default=TRUE,
              help="Run beta diversity envfit + Mantel [default %default]"),
  make_option(c("--run_maaslin3"), type="logical", default=TRUE,
              help="Run MaAsLin3 abundance+prevalence associations [default %default]"),
  make_option(c("--maaslin3_use_endpoint_only"), type="logical", default=TRUE,
              help="Use endpoint samples only in MaAsLin3 [default %default]"),
  make_option(c("--maaslin3_normalization"), type="character", default="TSS",
              help="MaAsLin3 normalization: TSS|CLR|NONE [default %default]"),
  make_option(c("--maaslin3_transform"), type="character", default="LOG",
              help="MaAsLin3 transform: LOG|PLOG|NONE [default %default]"),
  make_option(c("--covariates"), type="character", default=NULL,
              help="Comma-separated covariate names to include as fixed effects (must be columns in metadata). Optional."),
  make_option(c("--strata_for_mantel"), type="character", default=NULL,
              help="Column name to use as strata in Mantel permutations (e.g., Cohort). Optional."),
  make_option(c("--alpha_metrics"), type="character", default=NULL,
              help="Comma-separated list of alpha metric column names to test. If NULL, auto-detect from data."),
  make_option(c("--genus_counts_are_relative"), type="logical", default=FALSE,
              help="If TRUE, genus tables are already relative abundances and will not be TSS-normalized [default %default]")
)
parser <- OptionParser(option_list = opt_list, description = "Correlate tumor burden to metagenomic data (alpha, beta, taxa-level)")
opt <- parse_args(parser)

# ---------- Checks ----------
if (is.null(opt$metadata_files)) die("At least one --metadata_files path is required")
if (is.null(opt$genus_files)) die("At least one --genus_files path is required")
if (!file.exists(opt$tumor_file)) die(paste("Tumor file not found:", opt$tumor_file))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "intermediate"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "endpoint"), showWarnings = FALSE, recursive = TRUE)

set.seed(opt$seed)

# ---------- Helpers ----------
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) {
    readr::read_tsv(path, col_types = cols(.default = guess()))
  } else if (ext == "csv") {
    readr::read_csv(path, col_types = cols(.default = guess()))
  } else if (ext %in% c("xlsx","xls")) {
    readxl::read_xlsx(path)
  } else {
    die(sprintf("Unsupported file extension for %s", path))
  }
}

std_names <- function(df) {
  names(df) <- names(df) %>% stringr::str_replace_all("\\s+", "_") %>% tolower()
  df
}

require_cols <- function(df, cols, context) {
  miss <- setdiff(tolower(cols), names(df))
  if (length(miss) > 0) die(sprintf("Missing required column(s) in %s: %s", context, paste(miss, collapse=", ")))
}

# Metadata normalizer:
normalize_metadata <- function(df) {
  df <- std_names(df)
  if ("sample" %in% names(df)) {
    df <- dplyr::rename(df, Sample = sample)
  } else if ("group" %in% names(df)) {
    df <- dplyr::rename(df, Sample = group)
  } else {
    die("Metadata requires a 'Sample' (or 'sample'/'group') column")
  }
  if ("mouseid" %in% names(df)) df <- dplyr::rename(df, MouseID = mouseid)
  if ("cohort" %in% names(df))  df <- dplyr::rename(df, Cohort = cohort)
  if ("timepoint" %in% names(df)) df <- dplyr::rename(df, Timepoint = timepoint)
  require_cols(df, c("Sample","MouseID"), "metadata")
  if (!"Cohort" %in% names(df)) warnf("Metadata missing 'Cohort' column; downstream stratification may be limited")
  if (!"Timepoint" %in% names(df)) warnf("Metadata missing 'Timepoint' column; endpoint detection uses this by regex")
  df %>% dplyr::mutate(Sample = as.character(Sample),
                       MouseID = as.character(MouseID),
                       Cohort = Cohort %||% NA_character_,
                       Timepoint = Timepoint %||% NA_character_)
}

# Alpha normalizer:
normalize_alpha <- function(df) {
  df <- std_names(df)
  if ("sample" %in% names(df)) {
    df <- dplyr::rename(df, Sample = sample)
  } else if ("group" %in% names(df)) {
    df <- dplyr::rename(df, Sample = group)
  } else if ("id" %in% names(df)) {
    df <- dplyr::rename(df, Sample = id)
  } else {
    first <- names(df)[1]
    warnf(sprintf("Alpha file: assuming first column '%s' is Sample", first))
    df <- dplyr::rename(df, Sample = !!first)
  }
  df %>% dplyr::mutate(Sample = as.character(Sample))
}

# Genus table normalizer: genera rows, sample columns -> samples x genera
normalize_genus <- function(df) {
  df <- std_names(df)
  feat_col <- dplyr::case_when(
    "genus" %in% names(df) ~ "genus",
    "feature" %in% names(df) ~ "feature",
    TRUE ~ names(df)[1]
  )
  if (!(feat_col %in% names(df))) die("Genus table must have a first column of taxa names (e.g., 'Genus' or 'Feature')")
  df_long <- df %>% tidyr::pivot_longer(-all_of(feat_col), names_to = "Sample", values_to = "abund")
  df_wide <- df_long %>% tidyr::pivot_wider(names_from = all_of(feat_col), values_from = abund, values_fill = 0)
  df_wide %>% dplyr::mutate(Sample = as.character(Sample))
}

# TSS-normalize to relative abundances
tss_normalize <- function(df_samples_by_features, already_relative = FALSE) {
  if (already_relative) return(df_samples_by_features)
  mat <- df_samples_by_features %>% column_to_rownames("Sample") %>% as.matrix()
  rs <- rowSums(mat, na.rm=TRUE); rs[rs == 0] <- 1
  rel <- sweep(mat, 1, rs, FUN = "/")
  as.data.frame(rel, check.names = FALSE) %>% rownames_to_column("Sample")
}

# Tumor processing
parse_tp_index <- function(tp) {
  s <- as.character(tp)
  idx <- suppressWarnings(as.integer(stringr::str_extract(s, "[0-9]+")))
  ifelse(is.na(idx), NA_integer_, idx)
}
normalize_tumor <- function(df) {
  df <- std_names(df)
  if ("mouseid" %in% names(df)) df <- dplyr::rename(df, MouseID = mouseid)
  if ("cohort" %in% names(df)) df <- dplyr::rename(df, Cohort = cohort)
  if ("tumortimepoint" %in% names(df)) df <- dplyr::rename(df, TumorTimepoint = tumortimepoint)
  if ("tumorvolume_mm3" %in% names(df)) df <- dplyr::rename(df, TumorVolume_mm3 = tumorvolume_mm3)
  require_cols(df, c("MouseID","TumorTimepoint","TumorVolume_mm3"), "tumor file")
  df %>% dplyr::mutate(
    MouseID = as.character(MouseID),
    Cohort = ifelse("Cohort" %in% names(df), as.character(Cohort), NA_character_),
    TumorTP_index = parse_tp_index(TumorTimepoint),
    TumorVolume_mm3 = as.numeric(TumorVolume_mm3)
  )
}
get_tumor_endpoint <- function(tumor_df) {
  tumor_df %>%
    dplyr::group_by(MouseID) %>%
    dplyr::arrange(dplyr::desc(TumorTP_index), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(MouseID, Cohort, TumorTimepoint, TumorTP_index, TumorVolume_mm3) %>%
    dplyr::rename(TumorTimepoint_endpoint = TumorTimepoint,
                  TumorTP_index_endpoint = TumorTP_index,
                  TumorVolume_mm3_endpoint = TumorVolume_mm3)
}
get_microbiome_endpoint_samples <- function(meta_df, endpoint_regex) {
  if (!"Timepoint" %in% names(meta_df)) {
    warnf("Metadata has no Timepoint; taking last sample per MouseID by observed order")
    return(meta_df %>% group_by(MouseID) %>% slice(n()) %>% ungroup())
  }
  ep <- meta_df %>%
    dplyr::filter(!is.na(Timepoint) & stringr::str_detect(Timepoint, endpoint_regex)) %>%
    dplyr::group_by(MouseID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  fallback <- meta_df %>%
    dplyr::anti_join(ep, by = c("Sample","MouseID")) %>%
    dplyr::group_by(MouseID) %>%
    dplyr::slice(n()) %>%
    dplyr::ungroup()
  dplyr::bind_rows(ep, fallback) %>% dplyr::distinct(MouseID, .keep_all = TRUE)
}
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# ---------- Ingest ----------
inform("Reading tumor file...")
tumor_raw <- read_any(opt$tumor_file) %>% normalize_tumor()

inform("Reading metadata files...")
meta_list <- lapply(opt$metadata_files, function(p) { inform(paste0("  - ", p)); read_any(p) %>% normalize_metadata() })
metadata <- dplyr::bind_rows(meta_list) %>% dplyr::distinct()

inform("Reading alpha files (if any)...")
alpha <- NULL
if (!is.null(opt$alpha_files)) {
  alpha_list <- lapply(opt$alpha_files, function(p) { inform(paste0("  - ", p)); read_any(p) %>% normalize_alpha() })
  alpha <- dplyr::bind_rows(alpha_list) %>% dplyr::distinct()
}

inform("Reading genus tables...")
genus_list <- lapply(opt$genus_files, function(p) { inform(paste0("  - ", p)); read_any(p) %>% normalize_genus() })
genus_all <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Sample"), genus_list) %>% dplyr::mutate(across(-Sample, ~replace_na(., 0)))
genus_rel <- tss_normalize(genus_all, already_relative = isTRUE(opt$genus_counts_are_relative))

# ---------- Write masters ----------
write_tsv(metadata, file.path(opt$outdir, "intermediate", "master_metadata.tsv"))
if (!is.null(alpha)) write_tsv(alpha, file.path(opt$outdir, "intermediate", "master_alpha.tsv"))
write_tsv(genus_rel, file.path(opt$outdir, "intermediate", "master_genus_rel.tsv"))

# ---------- Build endpoint dataset ----------
inform("Building endpoint dataset (tumor endpoint + endpoint stool)...")
tumor_ep <- get_tumor_endpoint(tumor_raw)
meta_ep  <- get_microbiome_endpoint_samples(metadata, opt$endpoint_regex)

endpoint <- meta_ep %>%
  dplyr::left_join(tumor_ep, by = c("MouseID", "Cohort")) %>%
  dplyr::filter(!is.na(TumorVolume_mm3_endpoint))

if (nrow(endpoint) == 0) die("No endpoint matches between metadata and tumor endpoint. Check MouseID/Cohort matching and --endpoint_regex.")

# ---------- ALPHA: Spearman + optional LMM ----------
if (isTRUE(opt$run_alpha) && !is.null(alpha)) {
  inform("Running alpha diversity correlations (Spearman)...")
  alpha_cols <- if (!is.null(opt$alpha_metrics)) {
    strsplit(opt$alpha_metrics, ",")[[1]] %>% trimws()
  } else {
    setdiff(names(alpha), c("Sample"))
  }
  alpha_ep <- endpoint %>% dplyr::select(Sample, MouseID, Cohort, TumorVolume_mm3_endpoint) %>%
    dplyr::left_join(alpha, by = "Sample")
  res <- purrr::map_dfr(alpha_cols, function(m) {
    if (!m %in% names(alpha_ep)) { warnf(paste("Missing alpha metric:", m)); return(NULL) }
    s <- safe_spearman(alpha_ep[[m]], alpha_ep$TumorVolume_mm3_endpoint)
    tibble::tibble(metric = m, rho = s[["rho"]], p = s[["p"]], n = s[["n"]])
  }) %>% dplyr::mutate(q = p.adjust(p, method = "BH"))
  readr::write_tsv(res, file.path(opt$outdir, "endpoint", "alpha_vs_tumor.tsv"))

  if ("Cohort" %in% names(alpha_ep)) {
    inform("Fitting mixed models (alpha ~ scale(tumor) + (1|Cohort))...")
    lmm_res <- purrr::map_dfr(alpha_cols, function(m) {
      if (!m %in% names(alpha_ep)) return(NULL)
      dat <- alpha_ep %>% dplyr::select(all_of(m), TumorVolume_mm3_endpoint, Cohort) %>% tidyr::drop_na()
      if (nrow(dat) < 5) return(NULL)
      f <- stats::as.formula(sprintf("%s ~ scale(TumorVolume_mm3_endpoint) + (1|Cohort)", m))
      fit <- try(lme4::lmer(f, data = dat, REML = FALSE), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      broom.mixed::tidy(fit, effects = "fixed") %>%
        dplyr::filter(term == "scale(TumorVolume_mm3_endpoint)") %>%
        dplyr::mutate(metric = m) %>%
        dplyr::select(metric, estimate, std.error, statistic, p.value)
    })
    if (nrow(lmm_res) > 0) readr::write_tsv(lmm_res, file.path(opt$outdir, "endpoint", "alpha_lmm_vs_tumor.tsv"))
  }
} else if (isTRUE(opt$run_alpha)) {
  warnf("No alpha files were provided; skipping alpha analysis")
}

# ---------- BETA: envfit + Mantel ----------
if (isTRUE(opt$run_beta)) {
  inform("Running beta diversity analyses (Bray-Curtis + envfit + Mantel)...")
  keep_samples <- intersect(endpoint$Sample, genus_rel$Sample)
  if (length(keep_samples) < 3) die("<3 endpoint samples have genus data; cannot run beta analyses")
  X <- genus_rel %>% dplyr::filter(Sample %in% keep_samples) %>% tibble::column_to_rownames("Sample")
  X <- X[, colSums(X, na.rm=TRUE) > 0, drop=FALSE]
  meta_sub <- endpoint %>% dplyr::filter(Sample %in% keep_samples) %>% tibble::column_to_rownames("Sample")
  if (nrow(meta_sub) != nrow(X)) die("Sample mismatch in beta analysis")

  bray <- vegan::vegdist(X, method = "bray")
  ord  <- vegan::capscale(X ~ 1, distance = "bray")
  ef   <- vegan::envfit(ord ~ TumorVolume_mm3_endpoint, data = meta_sub %>% tibble::rownames_to_column("Sample"), permutations = opt$permutations)
  ef_tab <- tibble::tibble(
    variable = "TumorVolume_mm3_endpoint",
    r = as.numeric(ef$vectors$r),
    r2 = as.numeric(ef$vectors$r)^2,
    p = ef$vectors$pvals
  )
  readr::write_tsv(ef_tab, file.path(opt$outdir, "endpoint", "envfit_tumor.tsv"))

  tdist <- dist(scale(meta_sub$TumorVolume_mm3_endpoint), method = "euclidean")
  if (!is.null(opt$strata_for_mantel) && opt$strata_for_mantel %in% colnames(meta_sub)) {
    strata_vec <- as.factor(meta_sub[[opt$strata_for_mantel]])
    man <- vegan::mantel(bray, tdist, method = "spearman", permutations = opt$permutations, strata = strata_vec)
  } else {
    man <- vegan::mantel(bray, tdist, method = "spearman", permutations = opt$permutations)
  }
  mantel_tab <- tibble::tibble(statistic_r = unname(man$statistic), p = man$signif, permutations = opt$permutations)
  readr::write_tsv(mantel_tab, file.path(opt$outdir, "endpoint", "mantel_tumor.tsv"))
}

# ---------- TAXA-LEVEL: MaAsLin3 ----------
if (isTRUE(opt$run_maaslin3)) {
  inform("Running MaAsLin3 associations (abundance + prevalence models)...")

  features <- genus_rel
  meta2 <- metadata %>% dplyr::left_join(get_tumor_endpoint(tumor_raw), by = c("MouseID","Cohort"))

  if (isTRUE(opt$maaslin3_use_endpoint_only)) {
    samples_keep <- unique(endpoint$Sample)
    features <- features %>% dplyr::filter(Sample %in% samples_keep)
    meta2 <- meta2 %>% dplyr::filter(Sample %in% samples_keep)
  }

  if (!"TumorVolume_mm3_endpoint" %in% names(meta2)) die("TumorVolume_mm3_endpoint missing in prepared metadata for MaAsLin3")
  keep <- is.finite(meta2$TumorVolume_mm3_endpoint)
  meta2 <- meta2[keep, , drop=FALSE]
  features <- features %>% dplyr::filter(Sample %in% meta2$Sample)

  rownames(features) <- features$Sample; features$Sample <- NULL
  rownames(meta2) <- meta2$Sample; meta2$Sample <- NULL

  fixed_effects <- c("TumorVolume_mm3_endpoint")
  if (!is.null(opt$covariates)) {
    covs <- strsplit(opt$covariates, ",")[[1]] %>% trimws()
    covs <- covs[covs %in% colnames(meta2)]
    if (length(covs) > 0) fixed_effects <- c(fixed_effects, covs)
  }

  random_effects <- character(0)
  if ("Cohort" %in% colnames(meta2)) random_effects <- c(random_effects, "Cohort")
  if ("MouseID" %in% colnames(meta2)) random_effects <- c(random_effects, "MouseID")

  out_dir <- file.path(opt$outdir, "maaslin3_tumor")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  maaslin3::maaslin3(
    input_data   = as.data.frame(features, check.names = FALSE),
    input_metadata = as.data.frame(meta2, check.names = FALSE),
    output       = out_dir,
    fixed_effects = fixed_effects,
    random_effects = if (length(random_effects)>0) random_effects else NULL,
    normalization= opt$maaslin3_normalization,
    transform    = opt$maaslin3_transform,
    correction   = "BH",
    standardize  = TRUE,
    plot_summary_plot = TRUE,
    plot_associations = TRUE,
    cores = 1,
    save_models = FALSE,
    verbosity = "INFO"
  )
}

inform("All requested analyses finished successfully.")
