# Load required package
library(maaslin3)
# 1. Read data tables
genes <- read.delim("study.tsv",
                    row.names   = 1,
                    quote       = "",       # keep quotes in KO descriptions
                    check.names = FALSE)    # preserve sample IDs verbatim

metadata <- read.delim("study.txt",
                       row.names   = 1,
                       check.names = FALSE)
# 2. Align metadata rows to gene‑table columns
metadata <- metadata[colnames(genes), , drop = FALSE]
# 3. Cast categorical variables to factors
factor_levels <- list(
  sex                    = c("F", "M"),
  age                    = c("Y", "N"),
  percussion_sensitivity = c("Y", "N"),
  sinus_tract            = c("Y", "N"),
  probing                = c("Y", "N"),
  acceptable_restoration = c("Y", "N"),
  large_lesion           = c("Y", "N")
)

metadata[names(factor_levels)] <-
  Map(function(x, lv) factor(x, levels = lv),
      metadata[names(factor_levels)],
      factor_levels)
# 4. Run MaAsLin3
fit_out <- Maaslin3(
  input_data     = genes,
  input_metadata = metadata,
  output         = "maaslin3",
  formula        = "~ sex + age + percussion_sensitivity + sinus_tract +
                     probing + acceptable_restoration + large_lesion",
  max_pngs       = 100,
  cores          = 64
)
