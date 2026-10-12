# =============================================================================
#  01_feature_summary.R
#  Figure 2a (step 1) : collapse the per-sequence RNA secondary-structure
#  predictions into per-genomic-position means.
# =============================================================================
#
#  Input   : data/raw/SARS2_Structures_Detailed.csv       (large, not distributed)
#            columns: Genomic_Position, Sample_Index, GC, Energy, Structure
#  Output  : data/derived/SARS2_summary_results.csv
#            columns: Genomic_Position, Mean_GC, Mean_Energy, Dot_Proportion,
#                     Sample_Count
#
#  The Dot_Proportion is the mean fraction of unpaired (".") nucleotides in the
#  predicted structure at that position, which is later turned into the averaged
#  paired nucleotide proportion (Proppaired = 1 - Dot_Proportion).
#
#  Run from the repository root:
#      Rscript Figure_2/scripts/01_feature_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source("config.R")

file_path   <- file.path(DIR_DATA, "SARS2_Structures_Detailed.csv")
output_path <- file.path(DIR_DATA_DRV, "SARS2_summary_results.csv")

if (!file.exists(file_path)) {
  stop("Missing ", file_path,
       "\nThis is a ~300 MB per-sequence structure table; see data/README.md.")
}

# -- 1. memory-efficient read ------------------------------------------------
dt <- fread(file_path,
            select = c("Genomic_Position", "Sample_Index", "GC", "Energy", "Structure"),
            colClasses = list(character = "Structure"))

# -- 2. fraction of unpaired (".") positions ---------------------------------
calc_dot_prop <- function(x) {
  lengths <- nchar(x)
  dots <- sapply(regmatches(x, gregexpr("\\.", x)), length)
  mean(dots / lengths)
}

# -- 3. per-position summary -------------------------------------------------
result <- dt[, .(
  Mean_GC        = mean(GC, na.rm = TRUE),
  Mean_Energy    = mean(Energy, na.rm = TRUE),
  Dot_Proportion = calc_dot_prop(Structure),
  Sample_Count   = uniqueN(Sample_Index)
), by = Genomic_Position]

setorder(result, Genomic_Position)
fwrite(result, output_path)
message("wrote ", output_path, " (", nrow(result), " positions)")

rm(dt); gc()
