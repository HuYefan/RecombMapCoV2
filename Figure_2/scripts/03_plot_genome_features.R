# =============================================================================
#  03_plot_genome_features.R
#  Figure 2 : genomic feature distribution across the SARS-CoV-2 genome and the
#  cold / hot region comparison.
# =============================================================================
#
#  Panel a - genome-wide distribution of the averaged sequence and structural
#            features mapped alongside the averaged uniformed recombination
#            intensity (bottom track).
#  Panel b - cold regions (bottom 10% of recombination intensity) vs. hot
#            regions (top 10%), Mann-Whitney tests.
#
#  Input   : data/derived/SARS2_Results_260811.xlsx   (distributed)
#              columns: Genomic_Position_start/end, Mean_GC, Mean_Energy,
#                       Paired_Proportion, Pi_Site, Recombination_Intensity,
#                       Sample_Count, TRS
#  Output  : figures/main/Figure2a_genome_features.pdf
#            figures/main/Figure2b_cold_vs_hot_boxplots.pdf
#            figures/main/Figure2_genome_features_combined.pdf
#            figures/main/Figure2b_cold_hot_statistics.csv
#
#  Run from the repository root:
#      Rscript Figure_2/scripts/03_plot_genome_features.R
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("config.R")

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

# -- 2. cold / warm / hot classification -------------------------------------
# The manuscript uses the bottom / top 10% of the uniformed recombination
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

message(sprintf("cut-off %.0f%%: cold <= %.4f, hot >= %.4f  (%d cold / %d warm / %d hot)",
                CUTOFF * 100, cold_thresh, hot_thresh,
                sum(df_clean$Region_Type == "Cold"),
                sum(df_clean$Region_Type == "Warm"),
                sum(df_clean$Region_Type == "Hot")))

# -- 3. feature panels -------------------------------------------------------
# Mean_GC, Mean_Energy, Paired_Proportion, Pi_Site and the uniformed
# recombination intensity are the tracks shown in the published Figure 2.
# TRS is available in the feature table and is included in the Cold/Hot
# statistics table, but it is not one of the published tracks; set
# INCLUDE_TRS <- TRUE to add it as an extra facet.
INCLUDE_TRS <- FALSE

features <- c("Mean_GC", "Mean_Energy", "Paired_Proportion", "Pi_Site", "TRS")
extremes <- df_clean %>% filter(Region_Type != "Warm") %>% droplevels()

stats_tbl <- do.call(rbind, lapply(features, function(feat) {
  test <- wilcox.test(extremes[[feat]] ~ extremes$Region_Type)
  data.frame(
    Feature    = feat,
    Mean_Cold  = mean(extremes[[feat]][extremes$Region_Type == "Cold"], na.rm = TRUE),
    Mean_Hot   = mean(extremes[[feat]][extremes$Region_Type == "Hot"], na.rm = TRUE),
    W          = as.numeric(test$statistic),
    P_value    = test$p.value,
    row.names  = NULL
  )
}))
stats_tbl$P_adj_BH <- p.adjust(stats_tbl$P_value, method = "BH")
print(stats_tbl)

write.csv(stats_tbl, file.path(DIR_FIG_MAIN, "Figure2b_cold_hot_statistics.csv"),
          row.names = FALSE)

# -- 4. panel a : genome-wide feature tracks ---------------------------------
track_levels <- c("Mean_GC", "Mean_Energy", "Paired_Proportion", "Pi_Site")
if (INCLUDE_TRS) track_levels <- append(track_levels, "TRS", after = 4)
feature_levels <- c(track_levels, "Recombination_Intensity")

# explicit y-limits per track (keeps every facet on a comparable scale; the
# values are those used for the published figure)
track_limits <- list(
  Mean_GC                 = c(0, 1),
  Mean_Energy             = c(-80, -20),
  Paired_Proportion       = c(0, 1),
  Pi_Site                 = c(0, 0.004),
  TRS                     = c(0, 1.02),
  Recombination_Intensity = c(0, 15)
)
ylimits <- data.frame(
  Feature = feature_levels,
  x       = 0,
  y_min   = vapply(feature_levels, function(f) track_limits[[f]][1], numeric(1)),
  y_max   = vapply(feature_levels, function(f) track_limits[[f]][2], numeric(1))
) %>%
  mutate(Feature = factor(Feature, levels = feature_levels))

df_plot <- df_clean %>%
  select(Genomic_Position_start, all_of(track_levels), Recombination_Intensity) %>%
  pivot_longer(cols = -Genomic_Position_start, names_to = "Feature", values_to = "Value") %>%
  mutate(Genomic_Position_start = Genomic_Position_start + 49,
         Feature = factor(Feature, levels = feature_levels))

region_summary <- df_clean %>%
  mutate(run_id = cumsum(Region_Type != lag(Region_Type, default = first(Region_Type)))) %>%
  filter(Region_Type != "Warm") %>%
  group_by(run_id, Region_Type) %>%
  summarise(Start_Pos = min(Genomic_Position_start - 49),
            End_Pos   = max(Genomic_Position_start + 50),
            .groups = "drop") %>%
  select(Region_Type, Start_Pos, End_Pos) %>%
  arrange(Start_Pos)

p_genome <- ggplot(df_plot, aes(x = Genomic_Position_start, y = Value)) +
  geom_rect(data = region_summary,
            aes(xmin = Start_Pos, xmax = End_Pos, ymin = -Inf, ymax = Inf,
                fill = Region_Type),
            inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("Cold" = "dodgerblue", "Hot" = "firebrick")) +
  geom_line(color = "steelblue") +
  geom_blank(data = ylimits, aes(x = x, y = y_min), inherit.aes = FALSE) +
  geom_blank(data = ylimits, aes(x = x, y = y_max), inherit.aes = FALSE) +
  facet_grid(Feature ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(0, 29903), breaks = c(seq(0, 5) * 5000, 29903)) +
  labs(x = "Genomic Position", y = "Features") +
  theme_ms() +
  theme(legend.position = "bottom", legend.title = element_blank())

# -- 5. panel b : boxplots ---------------------------------------------------
box_limits <- list(
  Mean_GC                 = c(0.25, 0.6),
  Mean_Energy             = c(-80, -20),
  Paired_Proportion       = c(0.4, 0.8),
  Pi_Site                 = c(0, 0.00018),
  TRS                     = c(0, 0.001),
  Recombination_Intensity = c(0, 6)
)
ylimits2 <- data.frame(
  Feature = feature_levels,
  x       = "Hot",
  y_min   = vapply(feature_levels, function(f) box_limits[[f]][1], numeric(1)),
  y_max   = vapply(feature_levels, function(f) box_limits[[f]][2], numeric(1))
) %>%
  mutate(Feature = factor(Feature, levels = feature_levels))

df_long <- extremes %>%
  select(Region_Type, all_of(track_levels), Recombination_Intensity) %>%
  pivot_longer(cols = -Region_Type, names_to = "Feature", values_to = "Value") %>%
  mutate(Feature = factor(Feature, levels = feature_levels))

p_boxplot <- ggplot(df_long, aes(x = Region_Type, y = Value, fill = Region_Type)) +
  geom_boxplot(alpha = 0.7, outliers = FALSE) +
  scale_fill_manual(values = c("Cold" = "dodgerblue", "Hot" = "firebrick")) +
  geom_blank(data = ylimits2, aes(x = x, y = y_min), inherit.aes = FALSE) +
  geom_blank(data = ylimits2, aes(x = x, y = y_max), inherit.aes = FALSE) +
  facet_grid(Feature ~ ., scales = "free_y") +
  labs(x = "", y = "Feature Value") +
  theme_ms() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"))

# -- 6. save -----------------------------------------------------------------
# Geometry follows the published Figure 2 panels: the tracks and the boxplot
# column sit side by side on a 160 x 150 mm page, with the two panels also
# written separately.
FIG2_W <- 160
FIG2_H <- 150

gsave_pdf <- function(plot, file, width, height) {
  ggsave(plot, filename = file.path(DIR_FIG_MAIN, file),
         device = cairo_pdf, width = width, height = height, units = "mm")
  message("wrote ", file.path(DIR_FIG_MAIN, file))
}

final_plot <- (p_genome + p_boxplot) + plot_layout(widths = c(5.5, 1))
print(final_plot)
gsave_pdf(final_plot, "Figure2_genome_features_combined.pdf", FIG2_W, FIG2_H)
gsave_pdf(p_genome,  "Figure2a_genome_features.pdf", FIG2_W, FIG2_H)
gsave_pdf(p_boxplot, "Figure2b_cold_vs_hot_boxplots.pdf", FIG2_W / 5.5, FIG2_H)
