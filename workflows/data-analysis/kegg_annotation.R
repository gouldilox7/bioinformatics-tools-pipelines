#!/usr/bin/env Rscript
# ------------------------------------------------------------
# Annotate two HUMAnN tables for the Zapata dataset:
#   • study.txt      (UniRef → KO → defs & pathways)
#   • study.txt     (mapXXXXX → pathway titles)
#
# Output:
#   • study.tsv
#   • study.tsv
#
# Caches:
#   study.tsv, study.tsv, study.tsv
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(KEGGREST)
})

# -----------------------------------------------------------------
# Load UniRef90 → KO mapping
# -----------------------------------------------------------------
load_uniref_mapping <- function(
  map_file = "uniref90_ko.map.gz",
  map_url  = "https://bitbucket.org/biobakery/humann/raw/release_v3.0.1/maps/uniref90_ko.map.gz"
) {
  if (!file.exists(map_file)) {
    message("Downloading UniRef90→KO mapping...")
    download.file(map_url, destfile = map_file, mode = "wb")
  }
  read_tsv(gzfile(map_file), col_names = c("feature","KO"), show_col_types = FALSE)
}

# -----------------------------------------------------------------
# Helper: fetch in batches, skip 404s (used for pathway titles)
# -----------------------------------------------------------------
batch_get <- function(ids, prefix = "") {
  split(ids, ceiling(seq_along(ids)/10)) %>%
    map(~ {
      compact(map(.x, function(id) {
        out <- tryCatch(keggGet(paste0(prefix, id))[[1]], error = function(e) NULL)
        out
      }))
    }) %>% flatten()
}

# -----------------------------------------------------------------
# Annotate gene families: UniRef → KO → definition → pathways → title
# -----------------------------------------------------------------
annotate_genefamilies <- function(infile, outfile) {
  cat(">>> Reading gene families:", infile, "\n")
  df <- read_tsv(infile, show_col_types = FALSE)
  names(df)[1] <- "feature"
  # strip any taxonomy suffix
  df <- df %>% mutate(feature = sub("/path/to/input_dir|.*$", "", feature))

  # join mapping
  mapping <- load_uniref_mapping()
  df2 <- df %>% left_join(mapping, by = "feature") %>% filter(!is.na(KO))

  kos <- unique(df2$KO)
  cat("Mapped", nrow(df2), "rows to", length(kos), "unique KOs\n")

  # ---- KO → definition cache (bulk) ----
  existing_desc <- if (file.exists("study.tsv")) {
    read_tsv("study.tsv", show_col_types = FALSE)
  } else tibble(id=character(), desc=character())
  miss_kos <- setdiff(kos, existing_desc$id)
  if (length(miss_kos)) {
    cat("   ↳ Bulk fetching all KO definitions via keggList()\n")
    defs <- keggList("ko")
    new_desc <- tibble(
      id   = sub("^ko:", "", names(defs)),
      desc = unname(defs)
    ) %>% filter(id %in% miss_kos)
    existing_desc <- bind_rows(existing_desc, new_desc) %>% distinct()
    write_tsv(existing_desc, "study.tsv")
  }

  # ---- KO → pathway links cache ----
  ko2pw <- if (file.exists("study.tsv")) {
    read_tsv("study.tsv", show_col_types = FALSE)
  } else tibble(KO=character(), pathway=character())
  miss_links <- setdiff(kos, ko2pw$KO)
  if (length(miss_links)) {
    cat("   ↳ Fetching", length(miss_links), "KO→pathway links\n")
    links <- keggLink("pathway", paste0("ko:", miss_links))
    new_links <- tibble(
      KO      = sub("^ko:",   "", names(links)),
      pathway = sub("^path:", "", links)
    )
    ko2pw <- bind_rows(ko2pw, new_links) %>% distinct()
    write_tsv(ko2pw, "study.tsv")
  }

  # ---- pathway → title cache ----
  path_ids <- unique(ko2pw$pathway)
  pw_titles <- if (file.exists("study.tsv")) {
    read_tsv("study.tsv", show_col_types = FALSE)
  } else tibble(pathway=character(), pathway_name=character())
  miss_pw <- setdiff(path_ids, pw_titles$pathway)
  if (length(miss_pw)) {
    cat("   ↳ Fetching", length(miss_pw), "pathway titles\n")
    entries <- batch_get(miss_pw, prefix="path:")
    new_pw <- map_dfr(entries, ~ tibble(
      pathway      = sub("^path:", "", .$ENTRY),
      pathway_name = .$NAME
    ))
    pw_titles <- bind_rows(pw_titles, new_pw) %>% distinct()
    write_tsv(pw_titles, "study.tsv")
  }

  # merge & write
  out <- df2 %>%
    left_join(existing_desc, by=c("KO"="id")) %>%
    left_join(ko2pw,        by="KO") %>%
    left_join(pw_titles,    by="pathway")
  write_tsv(out, outfile)
  cat(">>> Written gene families (rows:", nrow(out), ")\n\n")
}

# -----------------------------------------------------------------
# Annotate pathway abundances: mapXXXXX → title
# -----------------------------------------------------------------
annotate_pathabundance <- function(infile, outfile) {
  cat(">>> Reading pathway abundances:", infile, "\n")
  df <- read_tsv(infile, show_col_types = FALSE)
  names(df)[1] <- "path_id"
  df <- df %>% filter(str_detect(path_id, "^map\\d{5}$"))
  pids <- str_extract(df$path_id, "\\d{5}$")

  # pathway → title cache
  pw_titles <- if (file.exists("study.tsv")) {
    read_tsv("study.tsv", show_col_types = FALSE)
  } else tibble(pathway=character(), pathway_name=character())
  miss_pw <- setdiff(pids, pw_titles$pathway)
  if (length(miss_pw)) {
    cat("   ↳ Fetching", length(miss_pw), "pathway titles\n")
    entries <- batch_get(miss_pw, prefix="path:")
    new_pw <- map_dfr(entries, ~ tibble(
      pathway      = sub("^path:", "", .$ENTRY),
      pathway_name = .$NAME
    ))
    pw_titles <- bind_rows(pw_titles, new_pw) %>% distinct()
    write_tsv(pw_titles, "study.tsv")
  }

  out <- df %>%
    mutate(pathway = pids) %>%
    left_join(pw_titles, by="pathway")
  write_tsv(out, outfile)
  cat(">>> Written pathway abundances (rows:", nrow(out), ")\n")
}

# -----------------------------------------------------------------
# Run annotations
# -----------------------------------------------------------------
annotate_genefamilies("study.txt",  "study.tsv")
annotate_pathabundance("study.txt","study.tsv")
