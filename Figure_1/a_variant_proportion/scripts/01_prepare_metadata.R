# =============================================================================
#  01_prepare_metadata.R
#  Figure 1a - step 1/5 : slice the GISAID metadata dump down to the columns
#  needed for the prevalence analysis (date-filtered, location-truncated).
# =============================================================================
#
#  Input   : data/raw/hcov-19_metadata1.tsv   
#  Output  : data/raw/hcov-19_metadata_filtered.csv
#            data/raw/hcov-19_AA_Substitutions.csv
#
#  Both outputs are large and are git-ignored.  Only "Collection date" of the
#  form YYYY-MM-DD is kept and "Location" is truncated from
#  "Country/State/..." to "Country", which is what the regional prevalence
#  panels (World / Asia / Europe / North America) use.
#
#  Run from the repository root:
#      Rscript Figure_1/a_variant_proportion/scripts/01_prepare_metadata.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

source("config.R")   # walks up from the working directory to the repo root

# -- input / output ----------------------------------------------------------
file_path <- file.path(DIR_DATA_RAW, "hcov-19_metadata.tsv")
out_meta  <- file.path(DIR_DATA_RAW, "hcov-19_metadata_filtered.csv")
out_aa    <- file.path(DIR_DATA_RAW, "hcov-19_AA_Substitutions.csv")

if (!file.exists(file_path)) {
  stop("Metadata dump not found: ", file_path,
       "\nSee data/README.md for how to obtain it.")
}

# -- 1. metadata slice -------------------------------------------------------
dt <- fread(file_path,
            select = c("Accession ID", "Collection date", "Location", "Pango lineage"),
            nThread = 4)

# keep only complete collection dates
dt <- dt[str_detect(`Collection date`, "^\\d{4}-\\d{2}-\\d{2}$")]

# Country/State/... -> Country
dt[, Location := str_trim(str_remove(Location, "/.*"))]

print(head(dt))
fwrite(dt, out_meta)
message("wrote ", out_meta, " (", nrow(dt), " records)")

# -- 2. amino-acid substitution slice (used by Figure S1) -------------------
dt_aa <- fread(file_path,
               select = c("Accession ID", "AA Substitutions", "Pango lineage"),
               quote = "",
               nThread = 4)
print(head(dt_aa))
fwrite(dt_aa, out_aa)
message("wrote ", out_aa, " (", nrow(dt_aa), " records)")
