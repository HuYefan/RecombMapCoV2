# =============================================================================
#  02_feature_windows.R
#  Figure 2 / Figure S3 - step 2/3 : sliding-window feature table.
# =============================================================================
#
#  Input   : data/derived/SARS2_summary_results.csv   (per-position means,
#            produced by 01_feature_summary.R)
#  Output  : data/derived/SARS2_Feature_Windows.csv
#
#  Window: 100 nt with a 10-nt step, the grid on which every feature is
#  averaged so that the sequence/structural tracks, nucleotide diversity and the
#  uniformed recombination intensity are directly comparable (Figure 2a,
#  Figure S3).
#
#  Running 01 + 02 rebuilds the feature side of the figure from the raw folding
#  and diversity outputs.  data/derived/SARS2_Results_260811.xlsx — the table
#  the published Figure 2 and Figure S3 were actually drawn from, with the
#  bottom/top 10% region classification already applied — ships with the
#  repository, so 03_plot_genome_features.R and
#  01_vcm_correlation_matrix.R run without re-running these two steps.
#
#  Run from the repository root:
#      Rscript Figure_2/scripts/02_feature_windows.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("config.R")

WINDOW_SIZE <- 100
STEP_SIZE   <- 10

# -- sliding-window mean over a set of value columns -------------------------
process_sliding_window <- function(data, pos_col, val_cols) {
  max_pos <- max(data[[pos_col]], na.rm = TRUE)
  starts  <- seq(1, max_pos - WINDOW_SIZE + 1, by = STEP_SIZE)

  map_df(starts, function(s) {
    e <- s + WINDOW_SIZE - 1
    data %>%
      filter(.data[[pos_col]] >= s, .data[[pos_col]] <= e) %>%
      summarise(across(all_of(val_cols), \(x) mean(x, na.rm = TRUE))) %>%
      mutate(window_start = s, window_end = e, window_mid = (s + e) / 2)
  }) %>%
    select(window_start, window_end, window_mid, everything())
}

# -- 1. sequence / structural features ---------------------------------------
summary_file <- file.path(DIR_DATA_DRV, "SARS2_summary_results.csv")
if (!file.exists(summary_file)) {
  stop("Missing ", summary_file,
       "\nRun Figure_2/scripts/01_feature_summary.R first, or provide ",
       "data/derived/SARS2_summary_results.csv (see data/README.md).")
}
feat <- read.csv(summary_file)
feat$Paired_Proportion <- 1 - feat$Dot_Proportion

features_windowed <- process_sliding_window(
  feat,
  pos_col  = "Genomic_Position",
  val_cols = c("Mean_GC", "Mean_Energy", "Paired_Proportion", "Sample_Count")
) %>%
  rename(GC_Ratio = Mean_GC, Min_Free_Energy = Mean_Energy)

# -- 2. nucleotide diversity -------------------------------------------------
entropy_file <- file.path(DIR_DATA_DRV, "SARS2_Genome_Diversity_Entropy.csv")
if (!file.exists(entropy_file)) {
  stop("Missing ", entropy_file, " - see data/README.md")
}
df_entropy <- read.csv(entropy_file)
entropy_windowed <- process_sliding_window(
  df_entropy,
  pos_col  = "Genomic_Position",
  val_cols = c("Pi_Site", "Entropy_Smoothed")
)

# -- 3. uniformed recombination intensity ------------------------------------
# published per-base track, read through the shared helper in config.R
df_recomb <- load_recomb_intensity()
recomb_windowed <- process_sliding_window(df_recomb, "pos", "intensity") %>%
  rename(Recombination_Intensity = intensity)

# -- 4. join on the shared window grid ---------------------------------------
final_combined <- features_windowed %>%
  inner_join(entropy_windowed, by = c("window_start", "window_end", "window_mid")) %>%
  inner_join(recomb_windowed, by = c("window_start", "window_end", "window_mid")) %>%
  mutate(Genomic_Position_start = window_start,
         Genomic_Position_end   = window_end)

out_file <- file.path(DIR_DATA_DRV, "SARS2_Feature_Windows.csv")
write.csv(final_combined, out_file, row.names = FALSE)
message("wrote ", out_file, " (", nrow(final_combined), " windows)")

print(head(final_combined))
