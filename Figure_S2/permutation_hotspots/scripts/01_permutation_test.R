# =============================================================================
#  01_permutation_test.R
#  Figure S2h-j : Monte-Carlo permutation test for genomic regions with
#  significantly elevated recombination intensity.
# =============================================================================
#
#  Method (matches the manuscript's statistics section)
#  ---------------------------------------------------
#  For each target region, the site-specific intensities inside the region are
#  compared against a background distribution comprising all other regions.
#  A one-sided permutation test (perm::permTS, Monte-Carlo, N = 5,000 resamples)
#  asks whether the mean intensity inside the region significantly exceeds the
#  background.  Raw p-values are corrected with Benjamini-Hochberg to keep the
#  FDR below 0.05.
#
#  Input   : data/reference/recomb_intensity_per_base.csv
#            data/reference/SARS-CoV-2_ORF.xlsx
#  Output  : figures/supplementary/FigureS2h_permutation_zscores.pdf
#            figures/supplementary/FigureS2h_permutation_results.tsv
#
#  Run from the repository root:
#      Rscript Figure_S2/permutation_hotspots/scripts/01_permutation_test.R
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(perm)
})

source("config.R")

set.seed(RNG_SEED)

# -- 1. inputs ---------------------------------------------------------------
# published per-base uniformed intensity track, read through config.R
recomb_data <- load_recomb_intensity()
orf_data <- as.data.frame(read_excel(file.path(DIR_DATA_REF, "SARS-CoV-2_ORF.xlsx")))

genome_vec <- recomb_data$intensity

# -- 2. permutation test -----------------------------------------------------
run_perm_ts_comparison <- function(r_start, r_end, n_mc = N_PERM) {
  indices_in <- max(1, r_start):min(length(genome_vec), r_end)
  x_scores <- genome_vec[indices_in]
  y_scores <- genome_vec[-indices_in]

  ts_res <- permTS(x = x_scores, y = y_scores,
                   alternative = "greater",
                   method = "exact.mc",
                   control = permControl(nmc = n_mc, seed = RNG_SEED))

  mean_x <- mean(x_scores); sd_x <- sd(x_scores)
  mean_y <- mean(y_scores); sd_y <- sd(y_scores)

  list(p_val = ts_res$p.value,
       z = (mean_x - mean_y) / sd(y_scores),
       obs_mean = mean_x, obs_sd = sd_x,
       bg_mean = mean_y, bg_sd = sd_y)
}

results_list <- list()
message("running permTS comparisons (", N_PERM, " resamples per region)...")
for (i in seq_len(nrow(orf_data))) {
  row <- orf_data[i, ]
  message("  testing ", row$ORF, " ...")
  results_list[[row$ORF]] <- run_perm_ts_comparison(row$Start, row$End)
}

summary_results <- tibble(
  ORF       = names(results_list),
  P_Value   = map_dbl(results_list, ~ .x$p_val),
  Z_Score   = map_dbl(results_list, ~ .x$z),
  Intensity = map_dbl(results_list, ~ .x$obs_mean),
  SD        = map_dbl(results_list, ~ .x$obs_sd)
) %>%
  mutate(padj = p.adjust(P_Value, method = "BH")) %>%
  arrange(desc(Z_Score))

out_tsv <- file.path(DIR_FIG_SUPP, "FigureS2h_permutation_results.tsv")
write.table(summary_results, out_tsv, sep = "\t", row.names = FALSE)
message("wrote ", out_tsv)

# -- 3. Z-score panel --------------------------------------------------------
p_res <- ggplot(summary_results, aes(x = reorder(ORF, Z_Score), y = Z_Score)) +
  geom_bar(stat = "identity", aes(fill = P_Value < 0.05), color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("grey", "firebrick"), name = "Significant\n(p < 0.05)") +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50") +
  labs(y = "Z-Score (Standard Deviations from Null Mean)", x = "Region") +
  theme_ms() +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2))

print(p_res)
ggsave(p_res,
       filename = file.path(DIR_FIG_SUPP, "FigureS2h_permutation_zscores.pdf"),
       device = cairo_pdf, width = 90, height = 100, units = "mm")

# -- 4. per-region intensity bars -------------------------------------------
plot_data <- summary_results %>%
  mutate(Is_Significant = ifelse(padj < 0.05, "Yes", "No"),
         ORF = factor(ORF, levels = unique(orf_data$ORF)))

p_bars <- ggplot(plot_data, aes(x = ORF, y = Intensity)) +
  geom_errorbar(aes(ymin = Intensity, ymax = Intensity + SD),
                width = 0.5, color = "black", linewidth = 0.5) +
  geom_col(aes(fill = Is_Significant), color = "black", width = 0.7, linewidth = 0.5) +
  scale_fill_manual(values = c(No = "grey", Yes = "firebrick"),
                    name = "Significant\n(p < 0.05)") +
  labs(x = "Region", y = "Average Recombination Intensity") +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "inside", legend.position.inside = c(0.2, 0.85))

print(p_bars)
ggsave(p_bars,
       filename = file.path(DIR_FIG_SUPP, "FigureS2h_permutation_intensity_bars.pdf"),
       device = cairo_pdf, width = 150, height = 60, units = "mm")
