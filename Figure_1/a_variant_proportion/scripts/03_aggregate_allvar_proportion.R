# =============================================================================
#  03_aggregate_allvar_proportion.R
#  Figure 1a - step 3/5 : same aggregation as step 2, but for the full panel of
#  variants of interest (Alpha ... XFG) used in the alternative prevalence
#  panel (all variants instead of recombinants only).
# =============================================================================
#
#  Input   : data/raw/gisaid_hcov-19_metadata_filtered_260401.csv
#            data/reference/Variants_of_Interest.xlsx
#  Output  : data/raw/gisaid_hcov-19_VOIproportion_260401.csv
#
#  Run from the repository root:
#      Rscript Figure_1/a_variant_proportion/scripts/03_aggregate_allvar_proportion.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(readxl)
})

source("config.R")

gisaid_file  <- file.path(DIR_DATA_RAW, "hcov-19_metadata_filtered.csv")
mapping_file <- file.path(DIR_DATA_REF, "Variants_of_Interest.xlsx")
output_file  <- file.path(DIR_DATA_RAW, "hcov-19_VOIproportion.csv")

if (!file.exists(gisaid_file)) {
  stop("Filtered metadata not found: ", gisaid_file,
       "\nRun 01_prepare_metadata.R first (see data/README.md).")
}

# -- 1. read the columns we need --------------------------------------------
gisaid_dt <- fread(gisaid_file,
                   sep = ",",
                   select = c("Collection date", "Location", "Pango lineage"),
                   check.names = FALSE)

# -- 2. expand the mapping (may or may not carry a Sublineage column) -------
mapping_raw <- as.data.table(read_excel(mapping_file))
if ("Sublineage" %in% names(mapping_raw)) {
  mapping_lookup <- mapping_raw[, list(
    Prefix = unlist(tstrsplit(Sublineage, "|", fixed = TRUE))
  ), by = list(Lineage)]
} else {
  mapping_lookup <- unique(mapping_raw[, list(Lineage, Prefix = Lineage)])
}
mapping_lookup[, Prefix := trimws(Prefix)]

# -- 3. prefix extraction ----------------------------------------------------
gisaid_dt[, Clean_Pango := sub(" \\(consensus call\\)", "", `Pango lineage`)]
gisaid_dt[, Prefix := sub("\\..*", "", Clean_Pango)]

# -- 4. join and label the remainder "Others" -------------------------------
results_dt <- merge(gisaid_dt, unique(mapping_lookup), by = "Prefix", all.x = TRUE)
results_dt[is.na(Lineage), Lineage := "Others"]

# -- 5. per date / location and global aggregates ---------------------------
stats_region <- results_dt[, list(count = .N),
                           by = list(`Collection date`, Location, Lineage)]
stats_region[, Total := sum(count), by = list(`Collection date`, Location)]
stats_region[, Proportion := count / Total]

stats_world <- results_dt[, list(count = .N), by = list(`Collection date`, Lineage)]
stats_world[, Location := "World"]
stats_world[, Total := sum(count), by = list(`Collection date`)]
stats_world[, Proportion := count / Total]

final_output <- rbind(stats_region, stats_world, use.names = TRUE)
setorder(final_output, `Collection date`, Location, -Proportion)

fwrite(final_output, output_file)
message("wrote ", output_file, " (", nrow(final_output), " rows)")

rm(gisaid_dt, results_dt); gc()
