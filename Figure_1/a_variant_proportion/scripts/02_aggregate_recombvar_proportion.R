# =============================================================================
#  02_aggregate_recombvar_proportion.R
#  Figure 1a - step 2/5 : aggregate variant counts and proportions per
#  collection date and location, for the recombinant-lineage panel.
# =============================================================================
#
#  Input   : data/raw/gisaid_hcov-19_metadata_filtered_260401.csv
#            data/reference/RecombVar_ofInterest.xlsx
#  Output  : data/raw/gisaid_hcov-19_proportion_260401.csv
#
#  Lineages are collapsed to their parent recombinant group through the
#  "Sublineage" mapping in RecombVar_ofInterest.xlsx (e.g. "XDV|NB" expands to
#  the prefixes XDV and NB).  Everything unmatched becomes "Others".  A global
#  "World" row set is appended next to the per-country rows.
#
#  Run from the repository root:
#      Rscript Figure_1/a_variant_proportion/scripts/02_aggregate_recombvar_proportion.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(readxl)
})

source("config.R")

gisaid_file <- file.path(DIR_DATA_RAW, "hcov-19_metadata_filtered.csv")
mapping_file <- file.path(DIR_DATA_REF, "RecombVar_ofInterest.xlsx")
output_file <- file.path(DIR_DATA_RAW, "hcov-19_proportion.csv")

if (!file.exists(gisaid_file)) {
  stop("Filtered metadata not found: ", gisaid_file,
       "\nRun 01_prepare_metadata.R first (see data/README.md).")
}

# -- 1. memory-efficient read of the columns we need -------------------------
gisaid_dt <- fread(gisaid_file,
                   sep = ",",
                   select = c("Collection date", "Location", "Pango lineage"),
                   check.names = FALSE)

# -- 2. expand the recombinant-group mapping ---------------------------------
mapping_raw <- as.data.table(read_excel(mapping_file))

# "XDV|NB" -> one row per prefix
mapping_lookup <- mapping_raw[, list(
  Prefix = unlist(tstrsplit(Sublineage, "|", fixed = TRUE))
), by = list(Lineage)]
mapping_lookup[, Prefix := trimws(Prefix)]

# -- 3. clean lineage strings and extract the prefix -------------------------
gisaid_dt[, Clean_Pango := sub(" \\(consensus call\\)", "", `Pango lineage`)]
gisaid_dt[, Prefix := sub("\\..*", "", Clean_Pango)]

# -- 4. left join to the mapping --------------------------------------------
results_dt <- merge(gisaid_dt, unique(mapping_lookup), by = "Prefix", all.x = TRUE)
results_dt[is.na(Lineage), Lineage := "Others"]

# -- 5. counts and proportions per date / location ---------------------------
stats_region <- results_dt[, list(count = .N),
                           by = list(`Collection date`, Location, Lineage)]
stats_region[, Total := sum(count), by = list(`Collection date`, Location)]
stats_region[, Proportion := count / Total]

# -- 6. global (World) aggregate --------------------------------------------
stats_world <- results_dt[, list(count = .N), by = list(`Collection date`, Lineage)]
stats_world[, Location := "World"]
stats_world[, Total := sum(count), by = list(`Collection date`)]
stats_world[, Proportion := count / Total]

# -- 7. combine and order ----------------------------------------------------
final_output <- rbind(stats_region, stats_world, use.names = TRUE)
setorder(final_output, `Collection date`, Location, -Proportion)

fwrite(final_output, output_file)
message("wrote ", output_file, " (", nrow(final_output), " rows)")

rm(gisaid_dt, results_dt); gc()
