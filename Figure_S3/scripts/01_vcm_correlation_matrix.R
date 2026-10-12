# =============================================================================
#  01_vcm_correlation_matrix.R
#  Figure S3 - correlation matrix (VCM) of the SARS-CoV-2 recombination
#  variables across distinct intensity regions.
# =============================================================================
#
#  Pairwise correlation analysis comparing the averaged uniformed recombination
#  intensity with the key sequence and structural features:
#    - RNA fragment genomic position
#    - averaged GC ratio            (Mean_GC)
#    - averaged minimum free energy (Mean_Energy, dG)
#    - averaged paired nucleotide proportion (Paired_Proportion)
#    - averaged nucleotide diversity (Pi_Site)
#
#  Scatter plots and Spearman correlation coefficients are shown for every
#  pairwise comparison, stratified by recombination region type
#  (Cold / Warm / Hot).
#
#  Input   : data/derived/SARS2_Results_260811.xlsx   (distributed)
#  Output  : figures/supplementary/FigureS3_vcm_correlation_matrix.pdf
#            figures/supplementary/FigureS3_pairwise_correlations.csv
#
#  Run from the repository root:
#      Rscript Figure_S3/scripts/01_vcm_correlation_matrix.R
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(GGally)
  library(ggplot2)
})

source("config.R")

set.seed(RNG_SEED)

# -- 1. read and clean -------------------------------------------------------
results_file <- file.path(DIR_DATA_DRV, "SARS2_Results_260811.xlsx")
if (!file.exists(results_file)) {
  stop("Missing ", results_file, " - see data/README.md")
}

df <- read_excel(results_file, col_names = TRUE)

cols_to_use <- c("Genomic_Position_start",
                 "Mean_GC", "Mean_Energy", "Paired_Proportion", "Pi_Site",
                 "TRS", "Recombination_Intensity")
df[cols_to_use] <- lapply(df[cols_to_use], as.numeric)
df_clean <- na.omit(df[cols_to_use])

# -- 2. cold / warm / hot labels (bottom / top 10%) --------------------------
CUTOFF <- 0.10
hot_thresh  <- quantile(df_clean$Recombination_Intensity, 1 - CUTOFF, na.rm = TRUE)
cold_thresh <- quantile(df_clean$Recombination_Intensity, CUTOFF, na.rm = TRUE)

df_clean <- df_clean %>%
  mutate(Region_Type = case_when(
    Recombination_Intensity >= hot_thresh  ~ "Hot",
    Recombination_Intensity <= cold_thresh ~ "Cold",
    TRUE                                   ~ "Warm"
  ),
  Region_Type = factor(Region_Type, levels = c("Cold", "Warm", "Hot")))

message(sprintf("cut-off %.0f%% -> %d Cold / %d Warm / %d Hot windows",
                CUTOFF * 100,
                sum(df_clean$Region_Type == "Cold"),
                sum(df_clean$Region_Type == "Warm"),
                sum(df_clean$Region_Type == "Hot")))

# -- 3. correlation matrix ---------------------------------------------------
vcm_vars <- c("Genomic_Position_start", "Mean_GC", "Mean_Energy",
              "Paired_Proportion", "Pi_Site", "Recombination_Intensity")

pairplot <- ggpairs(
  df_clean[, c(vcm_vars, "Region_Type")],
  aes(color = Region_Type, alpha = 0.5),
  upper = list(continuous = wrap("cor", method = "spearman", size = 2)),
  diag  = list(continuous = wrap("densityDiag")),
  lower = list(continuous = wrap("points", size = 0.5, stroke = 0)),
  title = "Variables Correlation Matrix"
) +
  scale_color_manual(values = c(Cold = "dodgerblue", Warm = "grey", Hot = "firebrick")) +
  scale_fill_manual(values = c(Cold = "dodgerblue", Warm = "grey", Hot = "firebrick")) +
  theme_bw() +
  theme(title = element_text(size = 8),
        text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 6),
        panel.grid.minor = element_blank())

print(pairplot)

out_pdf <- file.path(DIR_FIG_SUPP, "FigureS3_vcm_correlation_matrix.pdf")
ggsave(pairplot, filename = out_pdf, device = cairo_pdf,
       width = 180, height = 180, units = "mm")
message("wrote ", out_pdf)

# -- 4. numeric correlation table (global and per region type) ---------------
cor_for <- function(d, label) {
  vars <- c("Genomic_Position_start", "Mean_GC", "Mean_Energy",
            "Paired_Proportion", "Pi_Site")
  do.call(rbind, lapply(vars, function(v) {
    if (length(unique(d[[v]])) < 3) {
      return(data.frame(Stratum = label, Variable = v, Rho = NA_real_,
                        P_value = NA_real_, N = nrow(d)))
    }
    ct <- suppressWarnings(cor.test(d[[v]], d$Recombination_Intensity,
                                    method = "spearman", exact = FALSE))
    data.frame(Stratum = label, Variable = v,
               Rho = unname(ct$estimate), P_value = ct$p.value, N = nrow(d))
  }))
}

cor_tbl <- rbind(
  cor_for(df_clean, "All"),
  cor_for(df_clean %>% filter(Region_Type == "Cold"), "Cold"),
  cor_for(df_clean %>% filter(Region_Type == "Warm"), "Warm"),
  cor_for(df_clean %>% filter(Region_Type == "Hot"),  "Hot")
) %>%
  group_by(Stratum) %>%
  mutate(P_adj_BH = p.adjust(P_value, method = "BH")) %>%
  ungroup()

print(cor_tbl)
out_csv <- file.path(DIR_FIG_SUPP, "FigureS3_pairwise_correlations.csv")
write.csv(cor_tbl, out_csv, row.names = FALSE)
message("wrote ", out_csv)
