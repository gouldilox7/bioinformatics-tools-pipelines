#!/usr/bin/env Rscript
# =============================================================================
# significancetesting.R  —  One-stop CLI for nonparametric significance testing
# =============================================================================
# This script inspects your metadata and input table to pick appropriate tests automatically,
# runs them (plus post-hoc where appropriate), and writes clean CSVs named: <input_base>_<test>_results.csv
#
# Supported main tests (auto-picked per design; can be forced via --force):
#   - Independent 2-level: Mann–Whitney U / Wilcoxon rank-sum (stats::wilcox.test)
#   - Independent >2-levels: Kruskal–Wallis (stats::kruskal.test) + Dunn/Conover post-hoc (PMCMRplus)
#   - Repeated (paired) 2-level: Wilcoxon signed-rank (stats::wilcox.test)
#   - Repeated (blocked) >2-levels: Friedman (stats::friedman.test) + Conover post-hoc (PMCMRplus)
#   - Optional monotonic trend (ordered groups, independent): Jonckheere–Terpstra (clinfun::jonckheere.test)
#
# Effect sizes (added where well-defined):
#   - 2-group independent: Vargha–Delaney A (effsize::VD.A), Cliff’s delta (effsize::cliff.delta)
#   - Kruskal–Wallis: epsilon-squared (rcompanion::epsilonSquared)
#   - Friedman: Kendall’s W (rcompanion::kendallW; wraps DescTools::KendallW)
#   - 2-group paired: r for signed-rank (rcompanion::wilcoxonPairedR)
#
# ----------------------------------------------------------------------
# CLI flags (detailed, with generic examples; placeholders are generic)
# ----------------------------------------------------------------------
# --input (-i)         Path to the data table.
#                      Alpha mode: first column "sample", remaining columns numeric metrics.
#                      Abundance mode: first column feature ID, subsequent columns sample IDs.
#                      Example: --input metrics.tsv
#
# --metadata (-m)      Path to metadata; must contain a 'sample' column (case-insensitive).
#                      Example: --metadata meta.tsv
#
# --mode               auto|alpha|abundance  [default: auto]
#                      auto = detect layout from the input file structure.
#                      Example: --mode auto
#
# --id_col             (Optional) Metadata column for subjects/blocks (paired/repeated).
#                      Example: --id_col SubjectID
#
# --factors            (Optional) Comma-separated metadata columns to test; if omitted,
#                      all non-numeric, non-'sample', non-id_col columns are tested.
#                      Example: --factors Group,Treatment,Site
#
# --factor             (Optional) Comma-separated metadata columns to coerce to factor
#                      (useful for numeric-looking codes).
#                      Example: --factor Day,Batch
#
# --levels             (Optional) Set level orders: "Factor=A|B|C;Other=T0|T1|T2".
#                      Enables better ordering and trend tests when --trend is used.
#
# --nested             Parent>Child rules OR 'auto' OR 'off'  [default: manual (when string provided)]
#                      auto = infer parent→child rules with guardrails; off = disable nesting.
#                      Example: --nested auto
#                               --nested "Site>Group,Group>Treatment"
#
# --nested_min_parent_levels   Min parent levels that must qualify to enable A>B. [default: 1]
# --nested_min_child_levels    Min child levels required WITHIN each qualifying parent. [default: 2]
# --nested_min_per_group       Min samples per child level within parent (defaults to --min_per_group).
# --nested_max_rules           Cap the number of auto-detected rules. [default: 20]
# --nested_strict_only         If TRUE, only accept strict nesting (each child level under exactly one parent). [default: FALSE]
# --nested_prefer_coarse_parent If TRUE, require nlevels(parent) <= nlevels(child). [default: TRUE]
#
# --trend              (Flag) Also run Jonckheere–Terpstra trend tests for ordered factors (>=3 levels).
#                      Example: --trend
#
# --force              Force test: mannwhitney|ranksum|wilcoxon|kruskal|friedman|jt
#                      Example: --force kruskal
#
# --posthoc            Post-hoc for >2 levels: dunn|conover  [default: dunn]  (Friedman uses conover)
# --posthoc_alpha      Run post-hoc only if omnibus p <= alpha (set 1 to always). [default: 0.05]
#
# --p_adjust           p.adjust method for global FDR & PMCMRplus post-hoc: BH|holm|bonferroni|BY|... [default: BH]
# --fdr_scope          FDR scope for MAIN results: global|by_test [default: global]
#
# --subset             R expression on metadata to subset samples (NAs treated as FALSE).
#                      Example: --subset 'Cohort=="A" & Age>=18'
#
# --min_per_group      Minimum non-NA per group to run a test. [default: 2]
# --threads            Worker processes (future.apply). 1 = sequential. [default: 1]
# --diag               (Flag) Write diagnostic TSVs (counts per level per test; auto-nested proposals).
#
# OUTPUTS
# -------
#   <input_base>_mannwhitney_results.csv
#   <input_base>_wilcoxon_paired_results.csv
#   <input_base>_kruskal_results.csv
#   <input_base>_kruskal_posthoc_<dunn|conover>.csv
#   <input_base>_friedman_results.csv
#   <input_base>_friedman_posthoc_conover.csv
#   <input_base>_jt_results.csv
#
# ----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(rlang)
})

# ---- Robust dependency check helper -------------------------------------
ensure_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s. Install with install.packages('%s')",
                 pkg, why, pkg), call. = FALSE)
  }
}

# Required at runtime for core features:
ensure_pkg("matrixStats", "fast row/column summaries in abundance mode")
ensure_pkg("effsize", "effect sizes (Vargha–Delaney A, Cliff's delta)")
ensure_pkg("rcompanion", "epsilonSquared, Kendall's W, paired Wilcoxon effect size")
ensure_pkg("PMCMRplus", "Dunn/Conover & Friedman post-hoc tests")
ensure_pkg("clinfun", "Jonckheere–Terpstra trend test")

# ---- CLI ----------------------------------------------------------------
option_list <- list(
  make_option(c("-i","--input"), type="character", help="Path to input data table (required)."),
  make_option(c("-m","--metadata"), type="character", help="Path to metadata with a 'sample' column (required)."),
  make_option(c("--mode"), type="character", default="auto",
              help="auto|alpha|abundance [default: %default]"),
  make_option(c("--id_col"), type="character", default=NULL,
              help="Metadata column for subject/blocks (paired/repeated)."),
  make_option(c("--factors"), type="character", default=NULL,
              help="Comma-separated metadata columns to test (optional)."),
  make_option(c("--factor"), type="character", default=NULL,
              help="Comma-separated metadata columns to coerce to factor (optional)."),
  make_option(c("--levels"), type="character", default=NULL,
              help="Level order map, e.g. 'Group=Control|DrugA|DrugB;Treatment=Sham|Dose'."),
  make_option(c("--nested"), type="character", default=NULL,
              help="Parent>Child rules, OR 'auto', OR 'off'."),
  make_option(c("--nested_min_parent_levels"), type="integer", default=1,
              help="Min parent levels that qualify to enable A>B [default: %default]"),
  make_option(c("--nested_min_child_levels"), type="integer", default=2,
              help="Min child levels required WITHIN a parent level [default: %default]"),
  make_option(c("--nested_min_per_group"), type="integer", default=NA,
              help="Min samples per child level within parent (defaults to --min_per_group)"),
  make_option(c("--nested_max_rules"), type="integer", default=20,
              help="Max auto-detected nested rules to include [default: %default]"),
  make_option(c("--nested_strict_only"), action="store_true", default=FALSE,
              help="Only accept strict nesting (each child level under exactly one parent)"),
  make_option(c("--nested_prefer_coarse_parent"), type="logical", default=TRUE,
              help="Prefer coarser parents: require nlevels(parent) <= nlevels(child) [default: TRUE]"),
  make_option(c("--trend"), action="store_true", default=FALSE,
              help="Also run Jonckheere–Terpstra trend tests for ordered factors (>=3 levels)."),
  make_option(c("--force"), type="character", default=NULL,
              help="Force test: mannwhitney|ranksum|wilcoxon|kruskal|friedman|jt"),
  make_option(c("--posthoc"), type="character", default="dunn",
              help="Post-hoc (KW: dunn|conover; Friedman: conover) [default: %default]"),
  make_option(c("--posthoc_alpha"), type="double", default=0.05,
              help="Run post-hoc only if omnibus p <= this alpha (set 1 to always). [default: %default]"),
  make_option(c("--p_adjust"), type="character", default="BH",
              help="p.adjust method for global FDR and post-hoc (e.g. BH, holm, bonferroni, BY). [default: %default]"),
  make_option(c("--fdr_scope"), type="character", default="global",
              help="FDR scope for MAIN results: global|by_test [default: %default]"),
  make_option(c("--subset"), type="character", default=NULL,
              help="R expression on metadata to subset samples (e.g. 'Cohort==\"A\" & Age>=18')."),
  make_option(c("--min_per_group"), type="integer", default=2,
              help="Minimum non-NA per group to run a test [default: %default]"),
  make_option(c("--threads"), type="integer", default=1,
              help="Worker processes (future.apply). 1 = sequential [default: %default]"),
  make_option(c("--diag"), action="store_true", default=FALSE,
              help="Write diagnostic TSVs (counts per level per test; auto-nested proposals)."),
  make_option(c("--log"), type="character", default=NULL,
              help="Path to diagnostic log file (text). If omitted and --diag is set, defaults to 'auto_sigtests.log'.")
)
opt <- parse_args(OptionParser(option_list=option_list))
req <- c("input","metadata")
miss <- req[!req %in% names(opt) | sapply(opt[req], is.null)]
if (length(miss)) stop("Missing required arguments: ", paste(miss, collapse=", "))

# Parallel deps only if needed
if (opt$threads > 1) {
  ensure_pkg("future", "parallel execution with futures")
  ensure_pkg("future.apply", "parallel apply with futures")
}

# --- Diagnostic log setup ---
start_time <- Sys.time()
log_path <- opt$log
if (is.null(log_path) || !nzchar(log_path)) {
  if (isTRUE(opt$diag)) log_path <- "auto_sigtests.log"
}
log_line <- function(...) {
  if (!is.null(log_path) && nzchar(log_path)) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste(ts, "-", paste(..., collapse = " ")), "
", file = log_path, append = TRUE)
  }
}
pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error=function(e) NA_character_)
log_line("START auto_sigtests.R")
log_line("R:", R.version.string)
log_line("Packages:", paste0(
  "dplyr=", pkg_ver("dplyr"), ", tidyr=", pkg_ver("tidyr"), ", PMCMRplus=", pkg_ver("PMCMRplus"),
  ", effsize=", pkg_ver("effsize"), ", rcompanion=", pkg_ver("rcompanion"), ", clinfun=", pkg_ver("clinfun"),
  ", matrixStats=", pkg_ver("matrixStats"), ", future.apply=", pkg_ver("future.apply")))
log_line("Options parsed.")

# ---- Utilities -----------------------------------------------------------
basename_noext <- function(p) tools::file_path_sans_ext(basename(p))

parse_csv_or_tsv <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) readr::read_tsv(path, show_col_types = FALSE)
  else readr::read_csv(path, show_col_types = FALSE)
}

parse_levels_map <- function(s) {
  if (is.null(s) || s=="") return(list())
  parts <- str_split(s, ";")[[1]]
  out <- list()
  for (p in parts) {
    kv <- str_split_fixed(p, "=", 2)
    if (nchar(kv[1])==0 || nchar(kv[2])==0) next
    nm <- trimws(kv[1])
    levs <- str_split(kv[2], "\\|")[[1]] %>% trimws()
    out[[nm]] <- levs
  }
  out
}

eval_subset_mask <- function(expr_chr, data) {
  if (is.null(expr_chr) || expr_chr=="") return(rep(TRUE, nrow(data)))
  expr <- tryCatch(parse_expr(expr_chr), error = function(e) NULL)
  if (is.null(expr)) stop("Failed to parse --subset expression: ", expr_chr)
  out <- tryCatch(eval_tidy(expr, data = data), error = function(e) NULL)
  if (is.null(out) || !is.logical(out) || length(out)!=nrow(data))
    stop("--subset must evaluate to logical vector of nrow(metadata).")
  tidyr::replace_na(out, FALSE)
}

ordered_levels_from_data <- function(meta, col, mask, levels_map) {
  if (!is.null(levels_map[[col]])) return(levels_map[[col]])
  vals <- trimws(as.character(meta[[col]]))
  vals <- vals[mask & !is.na(vals)]
  unique(vals) # first appearance order
}

is_repeated_factor <- function(meta, id_col, f, mask) {
  if (is.null(id_col) || !id_col %in% names(meta)) return(FALSE)
  sub <- meta[mask & !is.na(meta[[f]]) & !is.na(meta[[id_col]]), c(id_col, f)]
  if (!nrow(sub)) return(FALSE)
  xt <- table(sub[[id_col]], sub[[f]])
  all(xt <= 1) && any(rowSums(xt > 0) >= 2)
}

# Effect size helpers
effsize_indep_2 <- function(x, g) {
  df <- tibble(x=x, g=g) %>% filter(!is.na(x) & !is.na(g))
  if (n_distinct(df$g)!=2 || nrow(df)<2) return(c(VD_A=NA_real_, CliffsDelta=NA_real_))
  vda <- tryCatch(effsize::VD.A(x ~ g, data=df)$estimate, error=function(e) NA_real_)
  cd  <- tryCatch(effsize::cliff.delta(x ~ g, data=df)$estimate, error=function(e) NA_real_)
  c(VD_A = as.numeric(vda), CliffsDelta = as.numeric(cd))
}

epsilon_sq_kw <- function(x, g) {
  df <- tibble(x=x, g=g) %>% filter(!is.na(x) & !is.na(g))
  if (n_distinct(df$g) < 2 || nrow(df)<2) return(NA_real_)
  tryCatch(as.numeric(rcompanion::epsilonSquared(df$x, df$g)), error=function(e) NA_real_)
}

kendall_w_friedman <- function(x, g, block, g_levels) {
  df <- tibble(x=x, g=factor(g, levels=g_levels), b=block) %>%
    filter(!is.na(x) & !is.na(g) & !is.na(b))
  ct <- df %>% count(b, g) %>% tidyr::pivot_wider(names_from=g, values_from=n, values_fill=0)
  complete_blocks <- ct$b[rowSums(as.matrix(ct[,-1])>0) == length(g_levels)]
  dfc <- df %>% filter(b %in% complete_blocks)
  if (!nrow(dfc)) return(NA_real_)
  wide <- tryCatch(tidyr::pivot_wider(dfc, names_from = b, values_from = x), error=function(e) NULL)
  if (is.null(wide)) return(NA_real_)
  mat <- as.matrix(wide[match(g_levels, wide$g), -1, drop=FALSE])
  tryCatch(as.numeric(rcompanion::kendallW(mat)), error=function(e) NA_real_)
}

add_global_fdr <- function(df, method, scope=c("global","by_test")) {
  scope <- match.arg(scope)
  if (!nrow(df) || !"P_value" %in% names(df)) return(df)
  if (scope=="by_test") {
    df %>% group_by(TestName) %>%
      mutate(FDR = p.adjust(P_value, method = method)) %>%
      ungroup()
  } else {
    df %>% mutate(FDR = p.adjust(P_value, method = method))
  }
}

# ---------- Auto-nesting inference ---------------------------------------
infer_nested_rules <- function(meta, factors, mask,
                               min_parent_levels = 1,
                               min_child_levels  = 2,
                               min_per_group     = 2,
                               strict_only       = FALSE,
                               prefer_coarse     = TRUE) {
  if (!length(factors)) return(tibble(parent=character(), child=character(),
                                      type=character(), n_parents_qual=integer()))
  rules <- list()
  meta2 <- meta[mask, , drop=FALSE]
  for (A in factors) {
    for (B in factors) {
      if (A == B) next
      pa <- trimws(as.character(meta2[[A]]))
      ch <- trimws(as.character(meta2[[B]]))
      ok <- !(is.na(pa) | is.na(ch))
      if (!any(ok)) next
      df <- tibble(A = pa[ok], B = ch[ok])
      if (!nrow(df)) next
      if (isTRUE(prefer_coarse)) {
        if (dplyr::n_distinct(df$A) > dplyr::n_distinct(df$B)) next
      }
      ct <- df %>% count(A, B, name = "n")
      ct2 <- ct %>% mutate(pass = n >= min_per_group) %>%
        group_by(A) %>% summarise(k_child = sum(pass), .groups="drop")
      qualifying_parents <- ct2 %>% filter(k_child >= min_child_levels) %>% pull(A)
      n_parents_qual <- length(qualifying_parents)
      if (n_parents_qual < min_parent_levels) next

      is_strict <- {
        mp <- ct %>% filter(n > 0) %>% group_by(B) %>% summarise(n_parents = n_distinct(A), .groups="drop")
        all(mp$n_parents == 1)
      }
      if (strict_only && !is_strict) next

      rules[[length(rules)+1]] <- tibble(parent = A, child = B,
                                         type = if (is_strict) "strict" else "stratified",
                                         n_parents_qual = n_parents_qual)
    }
  }
  if (!length(rules)) return(tibble(parent=character(), child=character(),
                                    type=character(), n_parents_qual=integer()))
  dplyr::bind_rows(rules) %>% dplyr::arrange(dplyr::desc(n_parents_qual), parent, child)
}

# ---- Load data & metadata -----------------------------------------------
dat  <- parse_csv_or_tsv(opt$input)
meta <- parse_csv_or_tsv(opt$metadata)

if (!"sample" %in% tolower(names(meta))) stop("Metadata must contain a 'sample' column (case-insensitive).")
names(meta)[tolower(names(meta))=="sample"] <- "sample"
meta$sample <- trimws(as.character(meta$sample))

if (!is.null(opt$factor) && nzchar(opt$factor)) {
  fc <- str_split(opt$factor, ",")[[1]] %>% trimws()
  for (cname in fc) if (cname %in% names(meta)) meta[[cname]] <- factor(meta[[cname]])
}

global_mask <- eval_subset_mask(opt$subset, meta)

# Determine factors to test (default: non-numeric columns except 'sample' and id_col)
if (is.null(opt$factors) || !nzchar(opt$factors)) {
  candidate <- setdiff(names(meta), c("sample", opt$id_col))
  fcols <- candidate[!sapply(meta[candidate], is.numeric)]
} else {
  fcols <- str_split(opt$factors, ",")[[1]] %>% trimws()
}
if (!all(fcols %in% names(meta))) {
  stop("Missing factor columns in metadata: ", paste(setdiff(fcols, names(meta)), collapse=", "))
}

levels_map <- parse_levels_map(opt$levels)
log_line("Factors:", paste(fcols, collapse=","))
log_line("id_col:", ifelse(length(opt$id_col) && !is.null(opt$id_col) && nzchar(opt$id_col), opt$id_col, "None"))
log_line("nested:", ifelse(is.null(opt$nested) || !nzchar(opt$nested), "manual/none", opt$nested))
log_line("trend:", as.character(opt$trend))
log_line("posthoc:", opt$posthoc, "alpha=", opt$posthoc_alpha, "p_adjust=", opt$p_adjust, "fdr_scope=", opt$fdr_scope)
log_line("min_per_group:", opt$min_per_group, "threads:", opt$threads)

# Validate p.adjust method
if (!opt$p_adjust %in% p.adjust.methods) {
  stop("--p_adjust must be one of: ", paste(p.adjust.methods, collapse=", "))
}

# ---- Detect mode & align -------------------------------------------------
mode <- tolower(opt$mode)
input_has_sample_col <- "sample" %in% tolower(names(dat))
samples_in_cols <- setdiff(names(dat), names(dat)[1])
metadata_samples <- meta$sample

if (mode == "auto") {
  if (input_has_sample_col) {
    mode <- "alpha"
  } else if (any(metadata_samples %in% samples_in_cols)) {
    mode <- "abundance"
  } else {
    stop("Could not auto-detect mode. Supply --mode alpha or --mode abundance.")
  }
} else if (!mode %in% c("alpha","abundance")) {
  stop("--mode must be auto|alpha|abundance")
}

if (mode == "alpha") {
  names(dat)[tolower(names(dat))=="sample"] <- "sample"
  if (!"sample" %in% names(dat)) stop("Alpha mode expects a 'sample' column in --input.")
  dat$sample <- trimws(as.character(dat$sample))
  keep <- dat %>% semi_join(meta, by="sample")
  meta2 <- meta %>% semi_join(keep, by="sample") %>% distinct(sample, .keep_all = TRUE)
  keep <- keep %>% inner_join(meta2 %>% select(sample), by="sample")
  metric_cols <- names(keep)[names(keep)!="sample"]
  keep <- keep %>% mutate(across(all_of(metric_cols), ~ suppressWarnings(as.numeric(.x))))
  A <- keep %>% arrange(match(sample, meta2$sample))
  meta_aligned <- meta2 %>% filter(sample %in% A$sample) %>% arrange(match(sample, A$sample))
  X_alpha <- as.matrix(A[, metric_cols, drop=FALSE])  # samples x metrics
  sample_vec <- A$sample
log_line("Mode:", mode, "| samples:", nrow(meta_aligned), "| metrics:", length(metric_cols))
} else {
  feat_ids <- dat[[1]]
  dat <- dat %>% mutate(across(-1, ~ suppressWarnings(as.numeric(.x))))
  all_samples <- colnames(dat)[-1]
  missing_in_input <- setdiff(metadata_samples, all_samples)
  if (length(missing_in_input)) {
    stop("These metadata samples are not columns in the feature table: ",
         paste(missing_in_input, collapse=", "))
  }
  meta2 <- meta %>% filter(sample %in% all_samples) %>% distinct(sample, .keep_all = TRUE)
  sample_order <- match(meta2$sample, all_samples)
  M <- as.matrix(dat[, -1, drop=FALSE])
  M <- M[, sample_order, drop=FALSE]
  feat_ids <- dat[[1]]
  meta_aligned <- meta2
}

# ---- Nested rules: manual | auto | off ----------------------------------
nested_mode <- if (is.null(opt$nested)) "manual" else {
  s <- tolower(trimws(opt$nested))
  if (s %in% c("auto","off")) s else "manual"
}

if (nested_mode == "off") {
  nested_rules <- tibble(parent=character(), child=character())
} else if (nested_mode == "auto") {
  npg <- if (is.na(opt$nested_min_per_group)) opt$min_per_group else opt$nested_min_per_group
  auto_rules <- infer_nested_rules(
    meta = meta_aligned,
    factors = fcols,
    mask = global_mask,
    min_parent_levels = opt$nested_min_parent_levels,
    min_child_levels  = opt$nested_min_child_levels,
    min_per_group     = npg,
    strict_only       = isTRUE(opt$nested_strict_only),
    prefer_coarse     = isTRUE(opt$nested_prefer_coarse_parent)
  )
  if (nrow(auto_rules) > opt$nested_max_rules) {
    auto_rules <- auto_rules %>% slice(seq_len(opt$nested_max_rules))
  }
  nested_rules <- auto_rules %>% select(parent, child)
  if (isTRUE(opt$diag)) {
    readr::write_tsv(auto_rules, "auto_sigtests_nested_rules.tsv")
  }
} else {
  parse_nested <- function(s) {
    if (is.null(s) || s=="") return(tibble(parent=character(), child=character()))
    pairs <- stringr::str_split(s, ",")[[1]] %>% trimws()
    pc <- stringr::str_split_fixed(pairs, ">", 2)
    tibble(parent=trimws(pc[,1]), child=trimws(pc[,2])) %>%
      dplyr::filter(parent!="" & child!="")
  }
  nested_rules <- parse_nested(opt$nested)
  if (nrow(nested_rules)) {
    if (!all(c(nested_rules$parent, nested_rules$child) %in% names(meta_aligned))) {
      stop("--nested references unknown metadata columns.")
    }
  }
}

log_line("Nested rules:", ifelse(exists("nested_rules") && nrow(nested_rules)>0, paste(nrow(nested_rules), "rules"), "none"))
# ---- Build Test Plan (Main + Nested) ------------------------------------
levels_for <- function(f) ordered_levels_from_data(meta_aligned, f, global_mask, levels_map)
tests <- list()

add_test <- function(name, factor_name, type, levs, mask, parent=NULL, parent_level=NULL) {
  tests[[length(tests)+1]] <<- list(
    name=name, factor=factor_name, type=type,
    levels=levs, mask=mask, parent=parent, parent_level=parent_level
  )
}

# Main effects
for (f in fcols) {
  levs <- levels_for(f)
  if (length(levs) < 2) next
  rep <- is_repeated_factor(meta_aligned, opt$id_col, f, global_mask)
  ttype <- if (!is.null(opt$force)) tolower(opt$force) else NULL
  if (!is.null(ttype)) {
    add_test(sprintf("%s (FORCED)", f), f, ttype, levs, global_mask)
    next
  }
  if (rep) {
    if (length(levs)==2) add_test(sprintf("%s: paired 2-level", f), f, "wilcoxon", levs, global_mask)
    if (length(levs)>=3) add_test(sprintf("%s: repeated >2-level", f), f, "friedman", levs, global_mask)
  } else {
    if (length(levs)==2) add_test(sprintf("%s: independent 2-level", f), f, "mannwhitney", levs, global_mask)
    if (length(levs)>=3) {
      add_test(sprintf("%s: independent >2-level", f), f, "kruskal", levs, global_mask)
      if (isTRUE(opt$trend)) add_test(sprintf("%s: monotonic trend", f), f, "jt", levs, global_mask)
    }
  }
}

# Nested: parent>child → run child within each parent level that exists
if (nrow(nested_rules)) {
  for (ri in seq_len(nrow(nested_rules))) {
    parent <- nested_rules$parent[ri]
    child  <- nested_rules$child[ri]
    p_levs <- levels_for(parent)
    for (pl in p_levs) {
      pmask <- global_mask & trimws(as.character(meta_aligned[[parent]])) == pl
      c_levs <- levels_for(child)
      rep <- is_repeated_factor(meta_aligned, opt$id_col, child, pmask)
      if (rep) {
        if (length(c_levs)==2) add_test(sprintf("%s|%s=%s", child, parent, pl), child, "wilcoxon", c_levs, pmask, parent, pl)
        if (length(c_levs)>=3) add_test(sprintf("%s|%s=%s", child, parent, pl), child, "friedman", c_levs, pmask, parent, pl)
      } else {
        if (length(c_levs)==2) add_test(sprintf("%s|%s=%s", child, parent, pl), child, "mannwhitney", c_levs, pmask, parent, pl)
        if (length(c_levs)>=3) {
          add_test(sprintf("%s|%s=%s", child, parent, pl), child, "kruskal", c_levs, pmask, parent, pl)
          if (isTRUE(opt$trend)) add_test(sprintf("%s trend|%s=%s", child, parent, pl), child, "jt", c_levs, pmask, parent, pl)
        }
      }
    }
  }
}

types <- sapply(tests, function(tt) tt$type)
type_tab <- as.data.frame(table(types), stringsAsFactors = FALSE)
log_line(paste0("Tests configured: ", length(tests), " | by type: ",
                paste(paste(type_tab$types, type_tab$Freq, sep="="), collapse=", ")))
if (!length(tests)) stop("No eligible tests were configured.")

# ---- Per-test runners ----------------------------------------------------
run_mannwhitney <- function(x, g, levs) {
  keep <- !is.na(x) & !is.na(g) & g %in% levs
  gk <- droplevels(factor(g[keep], levels=levs))
  xk <- x[keep]
  if (nlevels(gk)!=2) return(NULL)
  tab <- table(gk)
  if (any(tab < opt$min_per_group)) return(NULL)
  wt <- suppressWarnings(wilcox.test(x = xk[gk == levs[1]], y = xk[gk == levs[2]], paired = FALSE, exact = FALSE))
  meds <- tapply(xk, gk, median, na.rm=TRUE)
  eff  <- effsize_indep_2(xk, gk)
  tibble(
    Test="Mann-Whitney (Wilcoxon rank-sum)",
    DF=NA_integer_, Statistic=unname(wt$statistic), P_value=wt$p.value,
    N_A = as.integer(tab[levs[1]]), N_B = as.integer(tab[levs[2]]),
    Median_A = unname(meds[levs[1]]), Median_B = unname(meds[levs[2]]),
    VD_A = eff["VD_A"], CliffsDelta = eff["CliffsDelta"]
  )
}

run_wilcoxon_paired <- function(x, g, id, levs) {
  df <- tibble(x=x, g=g, id=id) %>% filter(!is.na(x) & !is.na(g) & !is.na(id) & g %in% levs)
  wide <- suppressWarnings(tidyr::pivot_wider(df, names_from=g, values_from=x))
  wide <- wide %>% select(id, all_of(levs))
  wide <- wide[stats::complete.cases(wide), , drop=FALSE]
  if (nrow(wide) < opt$min_per_group) return(NULL)
  wt <- suppressWarnings(wilcox.test(x = wide[[levs[1]]], y = wide[[levs[2]]],
                                   paired=TRUE, exact=FALSE))
  r_eff <- tryCatch({
    long_x <- as.numeric(t(as.matrix(wide[ , levs ])))
    long_g <- factor(rep(levs, each = nrow(wide)), levels = levs)
    as.numeric(rcompanion::wilcoxonPairedR(x = long_x, g = long_g))
  }, error=function(e) NA_real_)
  tibble(
    Test="Wilcoxon signed-rank (paired)",
    DF=NA_integer_, Statistic=unname(wt$statistic), P_value=wt$p.value,
    N_pairs = nrow(wide),
    Median_A = stats::median(wide[[levs[1]]], na.rm=TRUE),
    Median_B = stats::median(wide[[levs[2]]], na.rm=TRUE),
    Paired_r = r_eff
  )
}

run_kruskal <- function(x, g, levs) {
  keep <- !is.na(x) & !is.na(g) & g %in% levs
  gk <- droplevels(factor(g[keep], levels=levs))
  xk <- x[keep]
  if (nlevels(gk) < 2) return(NULL)
  tabs <- table(gk)
  if (any(tabs < opt$min_per_group)) return(NULL)
  kt <- suppressWarnings(kruskal.test(xk ~ gk))
  eps2 <- epsilon_sq_kw(xk, gk)
  tibble(
    Test="Kruskal-Wallis",
    DF=unname(kt$parameter), Statistic=as.numeric(kt$statistic), P_value=kt$p.value,
    N_groups = length(tabs), N_total = length(xk),
    EpsilonSquared = eps2
  )
}

run_friedman <- function(x, g, id, levs) {
  df <- tibble(x=x, g=factor(g, levels=levs), id=id) %>% filter(!is.na(x) & !is.na(g) & !is.na(id))
  ct <- df %>% count(id, g) %>% tidyr::pivot_wider(names_from=g, values_from=n, values_fill=0)
  ids_ok <- ct$id[rowSums(as.matrix(ct[,-1])>0) == length(levs)]
  dfc <- df %>% filter(id %in% ids_ok)
  if (!nrow(dfc)) return(NULL)
  ft <- suppressWarnings(friedman.test(x ~ g | id, data=dfc))
  W <- kendall_w_friedman(dfc$x, dfc$g, dfc$id, levs)
  tibble(
    Test="Friedman",
    DF=unname(ft$parameter), Statistic=as.numeric(ft$statistic), P_value=ft$p.value,
    N_blocks = n_distinct(dfc$id), N_groups = length(levs),
    KendallW = W
  )
}

run_jt <- function(x, g, levs) {
  keep <- !is.na(x) & !is.na(g) & g %in% levs
  gk <- droplevels(factor(g[keep], levels=levs, ordered=TRUE))
  xk <- x[keep]
  if (nlevels(gk) < 3) return(NULL)
  tabs <- table(gk)
  if (any(tabs < opt$min_per_group)) return(NULL)
  jt <- suppressWarnings(clinfun::jonckheere.test(xk, gk, alternative="two.sided"))
  tibble(
    Test="Jonckheere–Terpstra trend",
    DF=NA_integer_, Statistic=as.numeric(jt$statistic), P_value=jt$p.value,
    N_groups = length(tabs), N_total = length(xk)
  )
}

# Post-hoc helpers
posthoc_dunn <- function(x, g, levs, method) {
  df <- tibble(x=x, g=factor(g, levels=levs)) %>% filter(!is.na(x) & !is.na(g))
  if (n_distinct(df$g) < 2) return(NULL)
  ans <- suppressWarnings(PMCMRplus::kwAllPairsDunnTest(x ~ g, data=df, p.adjust.method = method))
  pv <- as.data.frame(ans$p.value); cn <- colnames(pv); rn <- rownames(pv)
  out <- list()
  for (i in seq_len(nrow(pv))) for (j in seq_len(ncol(pv))) if (!is.na(pv[i,j]) && i>j)
    out[[length(out)+1]] <- tibble(Level_A = rn[i], Level_B = cn[j], P_adj = as.numeric(pv[i,j]))
  if (!length(out)) return(NULL)
  bind_rows(out)
}

posthoc_conover_kw <- function(x, g, levs, method) {
  df <- tibble(x=x, g=factor(g, levels=levs)) %>% filter(!is.na(x) & !is.na(g))
  if (n_distinct(df$g) < 2) return(NULL)
  ans <- suppressWarnings(PMCMRplus::kwAllPairsConoverTest(x = df$x, g = df$g, p.adjust.method = method))
  pv <- as.data.frame(ans$p.value); cn <- colnames(pv); rn <- rownames(pv)
  out <- list()
  for (i in seq_len(nrow(pv))) for (j in seq_len(ncol(pv))) if (!is.na(pv[i,j]) && i>j)
    out[[length(out)+1]] <- tibble(Level_A = rn[i], Level_B = cn[j], P_adj = as.numeric(pv[i,j]))
  if (!length(out)) return(NULL)
  bind_rows(out)
}

posthoc_conover_friedman <- function(x, g, id, levs, method) {
  df <- tibble(x=x, g=factor(g, levels=levs), id=id) %>% filter(!is.na(x) & !is.na(g) & !is.na(id))
  ct <- df %>% count(id, g) %>% tidyr::pivot_wider(names_from=g, values_from=n, values_fill=0)
  ids_ok <- ct$id[rowSums(as.matrix(ct[,-1])>0) == length(levs)]
  dfc <- df %>% filter(id %in% ids_ok)
  if (!nrow(dfc)) return(NULL)
  ans <- suppressWarnings(PMCMRplus::frdAllPairsConoverTest(y = dfc$x, groups = dfc$g,
                                                            blocks = dfc$id, p.adjust.method = method))
  pv <- as.data.frame(ans$p.value); cn <- colnames(pv); rn <- rownames(pv)
  out <- list()
  for (i in seq_len(nrow(pv))) for (j in seq_len(ncol(pv))) if (!is.na(pv[i,j]) && i>j)
    out[[length(out)+1]] <- tibble(Level_A = rn[i], Level_B = cn[j], P_adj = as.numeric(pv[i,j]))
  if (!length(out)) return(NULL)
  bind_rows(out)
}

# ---- Run all tests across metrics/features -------------------------------
results_mw  <- list()
results_wxp <- list()
results_kw  <- list()
results_kw_ph <- list()
results_fr  <- list()
results_fr_ph <- list()
results_jt  <- list()

run_for_vector <- function(xvec, meta_rowmask, testdef) {
  f <- testdef$factor
  levs <- testdef$levels
  mask <- testdef$mask & meta_rowmask
  g <- trimws(as.character(meta_aligned[[f]]))
  g[!mask] <- NA_character_
  id <- if (!is.null(opt$id_col) && opt$id_col %in% names(meta_aligned)) meta_aligned[[opt$id_col]] else NA_character_

  out_main <- NULL; out_ph <- NULL

  if (testdef$type %in% c("mannwhitney","ranksum")) {
    out_main <- run_mannwhitney(xvec, g, levs)
    if (!is.null(out_main)) {
      out_main <- out_main %>% mutate(TestName = testdef$name, Factor=f,
                                      ParentFactor = testdef$parent %||% NA_character_,
                                      ParentLevel  = testdef$parent_level %||% NA_character_)
    }
  } else if (testdef$type=="wilcoxon") {
    out_main <- run_wilcoxon_paired(xvec, g, id, levs)
    if (!is.null(out_main)) {
      out_main <- out_main %>% mutate(TestName = testdef$name, Factor=f,
                                      ParentFactor = testdef$parent %||% NA_character_,
                                      ParentLevel  = testdef$parent_level %||% NA_character_)
    }
  } else if (testdef$type=="kruskal") {
    out_main <- run_kruskal(xvec, g, levs)
    if (!is.null(out_main)) {
      out_main <- out_main %>% mutate(TestName = testdef$name, Factor=f,
                                      ParentFactor = testdef$parent %||% NA_character_,
                                      ParentLevel  = testdef$parent_level %||% NA_character_)
      # NA-safe post-hoc trigger
      p <- out_main$P_value[1]
      if (isTRUE(opt$posthoc_alpha >= 1) || isTRUE(p <= opt$posthoc_alpha)) {
        if (tolower(opt$posthoc)=="dunn") out_ph <- posthoc_dunn(xvec, g, levs, opt$p_adjust)
        else                              out_ph <- posthoc_conover_kw(xvec, g, levs, opt$p_adjust)
        if (!is.null(out_ph)) out_ph <- out_ph %>%
            mutate(TestName = testdef$name, Factor=f,
                   ParentFactor = testdef$parent %||% NA_character_,
                   ParentLevel  = testdef$parent_level %||% NA_character_)
      }
    }
  } else if (testdef$type=="friedman") {
    out_main <- run_friedman(xvec, g, id, levs)
    if (!is.null(out_main)) {
      out_main <- out_main %>% mutate(TestName = testdef$name, Factor=f,
                                      ParentFactor = testdef$parent %||% NA_character_,
                                      ParentLevel  = testdef$parent_level %||% NA_character_)
      # NA-safe post-hoc trigger
      p <- out_main$P_value[1]
      if (isTRUE(opt$posthoc_alpha >= 1) || isTRUE(p <= opt$posthoc_alpha)) {
        out_ph <- posthoc_conover_friedman(xvec, g, id, levs, opt$p_adjust)
        if (!is.null(out_ph)) out_ph <- out_ph %>%
            mutate(TestName = testdef$name, Factor=f,
                   ParentFactor = testdef$parent %||% NA_character_,
                   ParentLevel  = testdef$parent_level %||% NA_character_)
      }
    }
  } else if (testdef$type=="jt") {
    out_main <- run_jt(xvec, g, levs)
    if (!is.null(out_main)) {
      out_main <- out_main %>% mutate(TestName = paste0(testdef$name," (JT)"), Factor=f,
                                      ParentFactor = testdef$parent %||% NA_character_,
                                      ParentLevel  = testdef$parent_level %||% NA_character_)
    }
  }

  list(main=out_main, posthoc=out_ph)
}

# Diagnostics
if (isTRUE(opt$diag)) {
  diag_rows <- list()
  for (tt in tests) {
    f <- tt$factor; levs <- tt$levels; mask <- tt$mask & global_mask
    g <- trimws(as.character(meta_aligned[[f]]))
    g[!mask] <- NA_character_
    tab <- table(factor(g, levels=levs), useNA="no")
    diag_rows[[length(diag_rows)+1]] <- tibble(
      TestName = tt$name, Factor = f, Type = tt$type,
      Level = names(tab), N = as.integer(tab)
    )
  }
  if (length(diag_rows)) bind_rows(diag_rows) %>% readr::write_tsv("auto_sigtests_diag_counts.tsv")
}

# Run per metric/feature
if (mode=="alpha") {
  metric_names <- colnames(X_alpha)
  meta_rowmask <- rep(TRUE, nrow(meta_aligned))
  runner <- function(j) {
    xvec <- X_alpha[, j]
    mains <- list(); phs <- list()
    for (tt in tests) {
      res <- run_for_vector(xvec, meta_rowmask, tt)
      if (!is.null(res$main)) mains[[length(mains)+1]] <- res$main %>% mutate(Feature = metric_names[j])
      if (!is.null(res$posthoc)) phs[[length(phs)+1]] <- res$posthoc %>% mutate(Feature = metric_names[j])
    }
    list(main = if (length(mains)) bind_rows(mains) else NULL,
         ph = if (length(phs)) bind_rows(phs) else NULL)
  }
  chunks <- if (opt$threads > 1) {
    future::plan(future::multisession, workers = opt$threads)
    on.exit(future::plan(future::sequential), add=TRUE)
    future.apply::future_lapply(seq_along(metric_names), runner)
  } else lapply(seq_along(metric_names), runner)

  for (chunk in chunks) {
    if (!is.null(chunk$main)) {
      df <- chunk$main
      for (i in seq_len(nrow(df))) {
        if (grepl("Mann-Whitney", df$Test[i])) results_mw[[length(results_mw)+1]] <- df[i,]
        else if (grepl("Wilcoxon signed", df$Test[i])) results_wxp[[length(results_wxp)+1]] <- df[i,]
        else if (grepl("Kruskal-Wallis", df$Test[i])) results_kw[[length(results_kw)+1]] <- df[i,]
        else if (grepl("Friedman", df$Test[i])) results_fr[[length(results_fr)+1]] <- df[i,]
        else if (grepl("Jonckheere", df$Test[i])) results_jt[[length(results_jt)+1]] <- df[i,]
      }
    }
    if (!is.null(chunk$ph)) {
      for (i in seq_len(nrow(chunk$ph))) {
        r <- chunk$ph[i,]
        results_kw_ph[[length(results_kw_ph)+1]] <- r
        results_fr_ph[[length(results_fr_ph)+1]] <- r
      }
    }
  }

} else {
  has_enough <- matrixStats::rowCounts(!is.na(M)) >= 2
  has_var    <- matrixStats::rowSds(M, na.rm = TRUE) > 0
  candidates <- which(has_enough & has_var)
log_line("Mode:", mode, "| samples:", ncol(M), "| features:", nrow(M), "| candidates:", length(candidates))
  runner <- function(i) {
    xvec <- M[i, ]
    mains <- list(); phs <- list()
    for (tt in tests) {
      res <- run_for_vector(xvec, rep(TRUE, length(xvec)), tt)
      if (!is.null(res$main)) mains[[length(mains)+1]] <- res$main %>% mutate(Feature = feat_ids[i])
      if (!is.null(res$posthoc)) phs[[length(phs)+1]] <- res$posthoc %>% mutate(Feature = feat_ids[i])
    }
    list(main = if (length(mains)) bind_rows(mains) else NULL,
         ph = if (length(phs)) bind_rows(phs) else NULL)
  }
  chunks <- if (opt$threads > 1) {
    future::plan(future::multisession, workers = opt$threads)
    on.exit(future::plan(future::sequential), add=TRUE)
    future.apply::future_lapply(candidates, runner)
  } else lapply(candidates, runner)

  for (chunk in chunks) {
    if (!is.null(chunk$main)) {
      df <- chunk$main
      for (i in seq_len(nrow(df))) {
        if (grepl("Mann-Whitney", df$Test[i])) results_mw[[length(results_mw)+1]] <- df[i,]
        else if (grepl("Wilcoxon signed", df$Test[i])) results_wxp[[length(results_wxp)+1]] <- df[i,]
        else if (grepl("Kruskal-Wallis", df$Test[i])) results_kw[[length(results_kw)+1]] <- df[i,]
        else if (grepl("Friedman", df$Test[i])) results_fr[[length(results_fr)+1]] <- df[i,]
        else if (grepl("Jonckheere", df$Test[i])) results_jt[[length(results_jt)+1]] <- df[i,]
      }
    }
    if (!is.null(chunk$ph)) {
      for (i in seq_len(nrow(chunk$ph))) {
        r <- chunk$ph[i,]
        results_kw_ph[[length(results_kw_ph)+1]] <- r
        results_fr_ph[[length(results_fr_ph)+1]] <- r
      }
    }
  }
}

# ---- Bind & write outputs ------------------------------------------------
input_base <- basename_noext(opt$input)

write_if_any <- function(df, fname) {
  if (!is.null(df) && nrow(df)) {
    readr::write_csv(df, fname, na = "")
    message("Wrote: ", fname, "  (", nrow(df), " rows)")
    log_line("Wrote:", fname, "rows:", nrow(df))
  }
}

finalize_main <- function(lst, family) {
  if (!length(lst)) return(NULL)
  df <- bind_rows(lst) %>%
    select(Feature, TestName, Factor, ParentFactor, ParentLevel,
           Test, DF, Statistic, P_value, everything())
  df <- add_global_fdr(df, method = opt$p_adjust, scope = tolower(opt$fdr_scope))
  df <- df %>% arrange(is.na(P_value), P_value)
  df
}

MW <- finalize_main(results_mw, "mannwhitney")
WXP <- finalize_main(results_wxp, "wilcoxon_paired")
KW <- finalize_main(results_kw, "kruskal")
FR <- finalize_main(results_fr, "friedman")
JT <- finalize_main(results_jt, "jt")

write_if_any(MW, sprintf("%s_mannwhitney_results.csv", input_base))
write_if_any(WXP, sprintf("%s_wilcoxon_paired_results.csv", input_base))
write_if_any(KW, sprintf("%s_kruskal_results.csv", input_base))
write_if_any(FR, sprintf("%s_friedman_results.csv", input_base))
write_if_any(JT, sprintf("%s_jt_results.csv", input_base))

if (!is.null(KW) && length(results_kw_ph)) {
  kw_names <- unique(KW$TestName)
  KWP <- bind_rows(results_kw_ph) %>% filter(TestName %in% kw_names) %>%
    select(Feature, TestName, Factor, ParentFactor, ParentLevel, Level_A, Level_B, P_adj) %>%
    arrange(is.na(P_adj), P_adj)
  write_if_any(KWP, sprintf("%s_kruskal_posthoc_%s.csv", input_base, tolower(opt$posthoc)))
}
if (!is.null(FR) && length(results_fr_ph)) {
  fr_names <- unique(FR$TestName)
  FRP <- bind_rows(results_fr_ph) %>% filter(TestName %in% fr_names) %>%
    select(Feature, TestName, Factor, ParentFactor, ParentLevel, Level_A, Level_B, P_adj) %>%
    arrange(is.na(P_adj), P_adj)
  write_if_any(FRP, sprintf("%s_friedman_posthoc_conover.csv", input_base))
}

if (!is.null(log_path) && nzchar(log_path)) log_line(sprintf("DONE in %.2f sec", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
message("✅ Done.")
