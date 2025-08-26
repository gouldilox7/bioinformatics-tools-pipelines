
# Load required libraries
library(tidyverse)
library(httr)
library(jsonlite)
library(parallel)
library(pbapply)

# ----- USER INPUT -----
input_files <- list.files(pattern = ".*\\.grouped\\.pathabundance\\.cpm\\.txt$")
output_dir <- "metacyc_api_outputs"
dir.create(output_dir, showWarnings = FALSE)

cache_file <- "study.tsv"
keywords <- c("toxin", "antibiotic", "resistance", "efflux", "multidrug")
batch_size <- 50  # number of pathway IDs per batch
num_cores <- 3    # conservative parallelism

# ----- Load or Initialize Cache -----
if (file.exists(cache_file)) {
  cache <- read.delim(cache_file, sep = "\t", stringsAsFactors = FALSE)
} else {
  cache <- tibble(Pathway = character(), Description = character())
}

# ----- BioCyc API Fetch Function -----
fetch_description_biocyc <- function(pid) {
  if (pid %in% cache$Pathway) {
    return(cache$Description[cache$Pathway == pid])
  }

  url <- paste0("/path/to/input_dir", pid, "/json")
  for (i in 1:3) {
    tryCatch({
      res <- GET(url)
      if (status_code(res) == 200) {
        content_json <- fromJSON(content(res, "text", encoding = "UTF-8"))
        name <- content_json$commonName
        return(name)
      }
    }, error = function(e) {})
    Sys.sleep(1)
  }
  return("Unknown")
}

# ----- Annotate Function for Each File -----
annotate_file <- function(file) {
  message("Processing file: ", file)
  df <- read.delim(file, sep = "\t", check.names = FALSE)
  colnames(df)[1] <- "Pathway"

  metacyc_df <- df %>%
    filter(grepl("^PWY|META|GLYCO|RXN", Pathway))

  pathway_ids <- unique(metacyc_df$Pathway)
  missing_ids <- setdiff(pathway_ids, cache$Pathway)

  if (length(missing_ids) > 0) {
    message("Missing IDs to fetch: ", length(missing_ids))
    batched_ids <- split(missing_ids, ceiling(seq_along(missing_ids) / batch_size))

    cl <- makeCluster(num_cores)
    for (i in seq_along(batched_ids)) {
      batch <- batched_ids[[i]]
      message("Fetching batch ", i, " of ", length(batched_ids), " (", length(batch), " IDs)...")
      clusterExport(cl, varlist = c("fetch_description_biocyc", "cache"), envir = environment())
      batch_results <- pbsapply(batch, fetch_description_biocyc, cl = cl)
      batch_cache <- tibble(Pathway = names(batch_results), Description = unname(batch_results))
      cache <<- bind_rows(cache, batch_cache)
      write.table(cache, cache_file, sep = "\t", row.names = FALSE, quote = FALSE)
      gc()
    }
    stopCluster(cl)
  }

  annotated_df <- left_join(metacyc_df, cache, by = "Pathway") %>%
    relocate(Description, .after = Pathway)

  filtered_df <- annotated_df %>%
    filter(str_detect(tolower(Description), paste(keywords, collapse = "|")))

  base <- tools::file_path_sans_ext(basename(file))
  write.table(annotated_df, file.path(output_dir, paste0(base, "study.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(filtered_df, file.path(output_dir, paste0(base, "study.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("✅ Finished: ", file)
}

# ----- Run for All Files -----
walk(input_files, annotate_file)
message("✅ All files processed. Outputs in: ", output_dir)
