# KW AUTO: main effects, optional nested within-tests, optional pairwise for >2 levels
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(stringr); library(tidyr); library(tibble)
  library(matrixStats); library(rlang)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Feature table (tab-delimited) [required]"),
  make_option(c("-m","--metadata"), type="character", help="Metadata with 'sample' + variables [required]"),
  make_option(c("-o","--output"), type="character", default="study.csv",
              help="Output CSV [default: %default]"),

  make_option(c("--factors"), type="character", help="Comma-separated metadata columns to test, e.g. 'Offspring,Treatment'"),
  make_option(c("--levels"), type="character", default=NULL,
              help="Optional factor level order, e.g. 'Offspring=Control|HF;Treatment=MDX|PHGG'"),
  make_option(c("--nested"), type="character", default=NULL,
              help="Optional nested rules, e.g. 'Offspring>Treatment,Center>Treatment'"),
  make_option(c("--subset"), type="character", default=NULL,
              help="Optional global subset evaluated on metadata, e.g. 'Site != \"Post /path/to/input_dir"'"),

  make_option(c("--pairwise"), type="logical", default=FALSE,
              help="If a factor has >2 levels, test all pairwise level pairs [default: %default]"),
  make_option(c("--min_per_group"), type="integer", default=2,
              help="Minimum non-NA per group to run a test [default: %default]"),
  make_option(c("--fdr_scope"), type="character", default="global",
              help="FDR scope: 'global' (all rows) or 'by_test' (per test name)"),
  make_option(c("--diag"), action="store_true", default=FALSE,
              help="Write diagnostics TSVs")
)
opt <- parse_args(OptionParser(option_list = option_list))
req <- c("input","metadata","factors")
miss <- req[!req %in% names(opt) | sapply(opt[req], is.null)]
if (length(miss)>0) stop("Missing required arguments: ", paste(miss, collapse=", "))

# ---------- Load & align ----------
feat <- read_tsv(opt$input, show_col_types = FALSE)
meta <- read_tsv(opt$metadata, show_col_types = FALSE)

# Normalize 'sample'
if (!"sample" %in% tolower(names(meta))) stop("Metadata must contain a 'sample' column (case-insensitive).")
names(meta)[tolower(names(meta))=="sample"] <- "sample"
meta$sample <- trimws(as.character(meta$sample))

# Clean feature column names & convert to numeric
colnames(feat) <- c(colnames(feat)[1], trimws(colnames(feat)[-1]))
feat <- feat %>% mutate(across(-1, ~ suppressWarnings(as.numeric(.x))))
all_samples <- colnames(feat)[-1]

missing_in_feature <- setdiff(meta$sample, all_samples)
if (length(missing_in_feature)>0) stop("These metadata samples are not columns in the feature table: ", paste(missing_in_feature, collapse=", "))

meta <- meta %>% filter(sample %in% all_samples) %>% distinct(sample, .keep_all = TRUE)
sample_order <- match(meta$sample, all_samples)
M <- as.matrix(feat[, -1, drop = FALSE])            # features x samples
M <- M[, sample_order, drop = FALSE]                 # reorder to metadata
feat_ids <- feat[[1]]

# ---------- Parse options ----------
factors <- str_split(opt$factors, ",", simplify = TRUE) %>% as.character() %>% trimws()
if (!all(factors %in% names(meta))) stop("Missing factor columns in metadata: ", paste(setdiff(factors, names(meta)), collapse=", "))

parse_levels_map <- function(s) {
  if (is.null(s) || is.na(s) || s == "") return(list())
  parts <- str_split(s, ";")[[1]]
  out <- vector("list", length(parts)); names(out) <- character(length(parts))
  j <- 0L
  for (p in parts) {
    kv <- str_split_fixed(p, "=", 2)
    if (nchar(kv[1]) == 0 || nchar(kv[2]) == 0) next
    j <- j + 1L
    nm <- trimws(kv[1])
    levs <- str_split(kv[2], "/path/to/input_dir|")[[1]] %>% trimws()
    out[[nm]] <- levs
    names(out)[j] <- nm
  }
  out[!sapply(out, is.null)]
}
levels_map <- parse_levels_map(opt$levels)

parse_nested <- function(s) {
  if (is.null(s) || is.na(s) || s == "") return(data.frame(parent=character(), child=character()))
  pairs <- str_split(s, ",")[[1]] %>% trimws()
  pc <- str_split_fixed(pairs, ">", 2)
  df <- tibble(parent = trimws(pc[,1]), child = trimws(pc[,2])) %>% filter(parent != "", child != "")
  df
}
nested_rules <- parse_nested(opt$nested)
if (nrow(nested_rules)) {
  if (!all(c(nested_rules$parent, nested_rules$child) %in% names(meta))) {
    stop("Nested rule references unknown metadata columns.")
  }
}

# Global subset mask
eval_subset <- function(expr_chr, data) {
  if (is.null(expr_chr) || is.na(expr_chr) || expr_chr == "") return(rep(TRUE, nrow(data)))
  expr <- tryCatch(parse_expr(expr_chr), error = function(e) NULL)
  if (is.null(expr)) stop("Failed to parse subset expression: ", expr_chr)
  out <- tryCatch(eval_tidy(expr, data = data), error = function(e) NULL)
  if (is.null(out) || !is.logical(out) || length(out) != nrow(data)) {
    stop("Subset did not evaluate to a logical vector of length nrow(metadata). Offending subset: ", expr_chr)
  }
  replace_na(out, FALSE)
}
global_mask <- eval_subset(opt$subset, meta)

# ---------- Helpers ----------
kw_two_group <- function(x, g, levA, levB, min_per_group=2) {
  keep <- !is.na(x) & !is.na(g) & (g == levA | g == levB)
  if (!any(keep)) return(NULL)
  gk <- droplevels(factor(g[keep], levels = c(levA, levB)))
  xk <- x[keep]
  if (nlevels(gk) < 2) return(NULL)
  tab <- table(gk)
  if (any(tab < min_per_group)) return(NULL)
  if (length(unique(xk)) < 2) return(NULL)
  kt <- tryCatch(kruskal.test(x = xk, g = gk), error = function(e) NULL)
  if (is.null(kt)) return(NULL)
  medA <- suppressWarnings(stats::median(xk[gk == levA], na.rm = TRUE))
  medB <- suppressWarnings(stats::median(xk[gk == levB], na.rm = TRUE))
  tibble(
    DF        = unname(kt$parameter),
    Statistic = as.numeric(kt$statistic),
    P_value   = kt$p.value,
    N_A = as.integer(tab[levA]), N_B = as.integer(tab[levB]),
    Median_A = medA, Median_B = medB
  )
}

# discover ordered levels from metadata (after global subset)
detect_levels <- function(col) {
  vals <- trimws(as.character(meta[[col]]))
  vals <- vals[global_mask & !is.na(vals)]
  unique(vals)  # keep first-appearance order
}

# Build test plan (list of tests to run)
tests <- list()
add_test <- function(name, factor_name, levA, levB, mask, test_type, parent=NULL, parent_level=NULL) {
  tests[[length(tests)+1L]] <<- list(
    name=name, factor=factor_name, A=levA, B=levB, mask=mask,
    type=test_type, parent=parent, parent_level=parent_level
  )
}

# 1) Main effects for each factor
for (f in factors) {
  levs <- if (!is.null(levels_map[[f]])) levels_map[[f]] else detect_levels(f)
  if (length(levs) < 2) next
  if (length(levs) == 2) {
    add_test(
      name = sprintf("%s: %s vs %s", f, levs[1], levs[2]),
      factor_name = f, levA = levs[1], levB = levs[2],
      mask = global_mask, test_type = "Main"
    )
  } else if (isTRUE(opt$pairwise)) {
    pr <- utils::combn(levs, 2, simplify = FALSE)
    for (p in pr) {
      add_test(
        name = sprintf("%s: %s vs %s", f, p[[1]], p[[2]]),
        factor_name = f, levA = p[[1]], levB = p[[2]],
        mask = global_mask, test_type = "Main"
      )
    }
  }
}

# 2) Nested rules: parent>child → run child comparisons within each parent level
if (nrow(nested_rules)) {
  for (ri in seq_len(nrow(nested_rules))) {
    parent <- nested_rules$parent[ri]; child <- nested_rules$child[ri]
    p_levs <- if (!is.null(levels_map[[parent]])) levels_map[[parent]] else detect_levels(parent)
    if (length(p_levs) < 1) next
    c_levs_all <- if (!is.null(levels_map[[child]])) levels_map[[child]] else detect_levels(child)

    for (pl in p_levs) {
      # mask: global subset AND parent == this level
      pmask <- global_mask & trimws(as.character(meta[[parent]])) == pl

      # detect child's levels within this parent subset if not provided
      c_levs <- if (is.null(levels_map[[child]])) {
        vals <- trimws(as.character(meta[[child]]))
        vals <- vals[pmask & !is.na(vals)]
        unique(vals)
      } else c_levs_all

      if (length(c_levs) < 2) next
      if (length(c_levs) == 2) {
        add_test(
          name = sprintf("%s: %s vs %s | %s=%s", child, c_levs[1], c_levs[2], parent, pl),
          factor_name = child, levA = c_levs[1], levB = c_levs[2],
          mask = pmask, test_type = "Nested", parent = parent, parent_level = pl
        )
      } else if (isTRUE(opt$pairwise)) {
        pr <- utils::combn(c_levs, 2, simplify = FALSE)
        for (p in pr) {
          add_test(
            name = sprintf("%s: %s vs %s | %s=%s", child, p[[1]], p[[2]], parent, pl),
            factor_name = child, levA = p[[1]], levB = p[[2]],
            mask = pmask, test_type = "Nested", parent = parent, parent_level = pl
          )
        }
      }
    }
  }
}

if (length(tests) == 0) {
  stop("No tests configured. Provide factors (and --pairwise if levels > 2), and/or nested rules.")
}

# ---------- Run ----------
# Pre-screen features: skip rows with no variance or too few observations overall
has_enough <- rowCounts(!is.na(M)) >= 2
has_var    <- rowSds(M, na.rm = TRUE) > 0
candidates <- which(has_enough & has_var)

# Prepare vectors used in tests
meta_vals <- lapply(factors, function(f) trimws(as.character(meta[[f]])))
names(meta_vals) <- factors

# We also need any child columns from nested that aren't in factors list
extra_children <- setdiff(unique(nested_rules$child), factors)
if (length(extra_children)) {
  extra_vals <- lapply(extra_children, function(f) trimws(as.character(meta[[f]])))
  names(extra_vals) <- extra_children
  meta_vals <- c(meta_vals, extra_vals)
}
# And any parents not in factors
extra_parents <- setdiff(unique(nested_rules$parent), names(meta_vals))
if (length(extra_parents)) {
  extra_vals2 <- lapply(extra_parents, function(f) trimws(as.character(meta[[f]])))
  names(extra_vals2) <- extra_parents
  meta_vals <- c(meta_vals, extra_vals2)
}

build_group <- function(test) {
  g <- trimws(as.character(meta[[test$factor]]))
  # apply mask by setting others to NA so they drop out
  g[!test$mask] <- NA_character_
  g
}

# Optionally write diagnostics of counts per test
if (isTRUE(opt$diag)) {
  diag_list <- vector("list", length(tests))
  for (ti in seq_along(tests)) {
    tt <- tests[[ti]]
    g <- build_group(tt)
    tab <- table(factor(g, levels = c(tt$A, tt$B)), useNA = "no")
    diag_list[[ti]] <- tibble(
      Test = tt$name, Factor = tt$factor, Type = tt$type,
      Level_A = tt$A, Level_B = tt$B,
      N_A = as.integer(tab[tt$A]), N_B = as.integer(tab[tt$B])
    )
  }
  bind_rows(diag_list) %>% write_tsv("study.tsv")
}

out <- vector("list", length(candidates) * length(tests))
k <- 0L

for (i in candidates) {
  x <- M[i, ]
  for (ti in seq_along(tests)) {
    tt <- tests[[ti]]
    g <- build_group(tt)
    res <- kw_two_group(x, g, tt$A, tt$B, min_per_group = opt$min_per_group)
    if (!is.null(res)) {
      k <- k + 1L
      out[[k]] <- tibble(
        Feature      = feat_ids[i],
        TestName     = tt$name,
        TestType     = tt$type,
        Factor       = tt$factor,
        ParentFactor = ifelse(is.null(tt$parent), NA_character_, tt$parent),
        ParentLevel  = ifelse(is.null(tt$parent_level), NA_character_, tt$parent_level),
        Level_A      = tt$A,
        Level_B      = tt$B,
        Test         = "Kruskal-Wallis"
      ) %>% bind_cols(res)
    }
  }
}

res <- bind_rows(out)
if (nrow(res) == 0) {
  message("No valid tests were run. Causes: too few per-group samples, missing levels, or no variance.")
  quit(status = 0)
}

# FDR
if (tolower(opt$fdr_scope) == "by_test") {
  res <- res %>%
    group_by(TestName) %>%
    mutate(FDR_adjusted_p = p.adjust(P_value, method = "BH")) %>%
    ungroup() %>%
    arrange(TestName, FDR_adjusted_p, P_value)
} else {
  res <- res %>%
    mutate(FDR_adjusted_p = p.adjust(P_value, method = "BH")) %>%
    arrange(FDR_adjusted_p, P_value)
}

res <- res %>%
  select(Feature, TestName, TestType, Factor, ParentFactor, ParentLevel,
         Level_A, Level_B, Test, DF, Statistic, P_value, FDR_adjusted_p,
         N_A, N_B, Median_A, Median_B)

write_csv(res, opt$output)
cat("✅ Wrote ", nrow(res), " rows to ", opt$output, "\n", sep = "")
