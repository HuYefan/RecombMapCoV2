# =============================================================================
#  04_plot_recombvar_proportion.R
#  Figure 1a - step 4/5 : stacked area charts of recombinant-lineage
#  prevalence (XBB, XEC, XFG, XDV vs. everything else) for World, Asia, Europe
#  and North America, using a 7-day centred rolling mean.
# =============================================================================
#
#  Input   : data/raw/gisaid_hcov-19_proportion_260401.csv   (distributed, 1.8 MB)
#  Output  : figures/main/Figure1a_recombvar_proportion_area.pdf
#            figures/main/Figure1a_recombvar_proportion_line.pdf
#
#  Run from the repository root:
#      Rscript Figure_1/a_variant_proportion/scripts/04_plot_recombvar_proportion.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggplot2)
  library(slider)   # rolling means
})

source("config.R")

# -- 1. data ----------------------------------------------------------------
data <- read_csv(file.path(DIR_DATA_RAW, "hcov-19_proportion.csv"),
                 show_col_types = FALSE)

DATE_MIN <- as.Date("2019-12-01")
DATE_MAX <- as.Date("2026-04-01")

plot_data <- data %>%
  mutate(Date = as.Date(`Collection date`)) %>%
  filter(Date >= DATE_MIN, Date <= DATE_MAX,
         Location %in% c("World", "Asia", "Europe", "North America")) %>%
  complete(Date = seq.Date(DATE_MIN, DATE_MAX, by = "day"),
           Location, Lineage,
           fill = list(count = 0, Proportion = 0)) %>%
  mutate(
    Proportion = replace(Proportion, Proportion < 0, 0),
    Variant = factor(Lineage,
                     levels = c("Others", "XBB", "XEC", "XFG", "XDV", "OtherRecombVar"),
                     labels = c("Others", "XBB*", "XEC*", "XFG*", "XDV*", "OtherRecombVar")),
    Region = factor(Location,
                    levels = c("World", "Asia", "Europe", "North America"))
  ) %>%
  arrange(Region, Variant, Date) %>%
  group_by(Region, Variant) %>%
  mutate(
    # 7-day centred rolling sum (matching the manuscript's smoothing)
    Smoothed_count = slide_index_dbl(count, Date, sum,
                                     .before = days(3), .after = days(3)),
    Smoothed_Total = slide_index_dbl(Total, Date, sum,
                                     .before = days(3), .after = days(3)),
    Smoothed_Proportion = Smoothed_count / Smoothed_Total
  ) %>%
  ungroup() %>%
  select(Date, Variant, Region, Proportion, Smoothed_Proportion)

# -- 2. colour map ----------------------------------------------------------
VariantList <- c("Others", "XBB*", "XEC*", "XFG*", "XDV*", "OtherRecombVar")
ColorList   <- c("grey90", "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "grey")
Color_mapping <- setNames(ColorList, VariantList)

x_breaks <- as.Date(c("2020-01-01", "2020-07-01",
                      "2021-01-01", "2021-07-01",
                      "2022-01-01", "2022-07-01",
                      "2023-01-01", "2023-07-01",
                      "2024-01-01", "2024-07-01",
                      "2025-01-01", "2025-07-01",
                      "2026-01-01"))

# -- 3. stacked area chart --------------------------------------------------
p1 <- ggplot(plot_data, aes(x = Date, y = Smoothed_Proportion, fill = Variant)) +
  geom_area(linewidth = 0.2, colour = "black") +
  facet_wrap(~ Region, ncol = 1) +
  scale_x_date(date_labels = "%Y-%m", breaks = x_breaks,
               limits = c(DATE_MIN, DATE_MAX), expand = expansion(mult = 0)) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = Color_mapping) +
  labs(x = "Date", y = "Proportion", fill = "Variant") +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = NA, color = NA)
  )

print(p1)
ggsave(p1, filename = file.path(DIR_FIG_MAIN, "Figure1a_recombvar_proportion_area.pdf"),
       device = cairo_pdf, width = 60, height = 120, units = "mm")

# -- 4. line version (recombinants only) ------------------------------------
p2 <- ggplot(plot_data %>% filter(Variant != "Others"),
             aes(x = Date, y = Smoothed_Proportion, color = Variant)) +
  geom_line(linewidth = 0.2) +
  facet_wrap(~ Region, ncol = 1) +
  scale_x_date(date_labels = "%Y-%m", breaks = x_breaks,
               limits = c(DATE_MIN, DATE_MAX), expand = expansion(mult = 0)) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = Color_mapping) +
  labs(x = "Date", y = "Proportion", color = "Variant") +
  theme_bw() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = NA, color = NA)
  )

print(p2)
ggsave(p2, filename = file.path(DIR_FIG_MAIN, "Figure1a_recombvar_proportion_line.pdf"),
       device = cairo_pdf, width = 60, height = 120, units = "mm")
