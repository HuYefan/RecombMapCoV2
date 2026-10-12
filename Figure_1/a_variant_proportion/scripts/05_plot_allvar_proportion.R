# =============================================================================
#  05_plot_allvar_proportion.R
#  Figure 1a - step 5/5 : stacked area chart of all variants of interest
#  (Alpha ... XFG) over time, World / Asia / Europe / North America.
# =============================================================================
#
#  Input   : data/raw/hcov-19_VOIproportion.csv   (distributed, 3.0 MB)
#  Output  : figures/main/Figure1a_allvar_proportion_area.pdf
#
#  Run from the repository root:
#      Rscript Figure_1/a_variant_proportion/scripts/05_plot_allvar_proportion.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggplot2)
  library(slider)
})

source("config.R")

# -- 1. data ----------------------------------------------------------------
data <- read_csv(file.path(DIR_DATA_RAW, "hcov-19_VOIproportion.csv"),
                 show_col_types = FALSE)

DATE_MIN <- as.Date("2019-12-01")
DATE_MAX <- as.Date("2026-03-31")

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
                     levels = c("Alpha", "Delta", "Omicron", "BA.1", "BA.2", "BA.3",
                                "BA.4", "BA.5", "XBB", "JN.1", "XEC", "XFG", "XDV",
                                "OtherRecombVar", "Others"),
                     labels = c("Alpha", "Delta", "Omicron", "BA.1", "BA.2", "BA.3",
                                "BA.4", "BA.5", "XBB", "JN.1", "XEC", "XFG", "XDV",
                                "OtherRecombVar", "OtherVar")),
    Region = factor(Location,
                    levels = c("World", "Asia", "Europe", "North America"))
  ) %>%
  arrange(Region, Variant, Date) %>%
  group_by(Region, Variant) %>%
  mutate(
    Smoothed_count = slide_index_dbl(count, Date, sum,
                                     .before = days(3), .after = days(3)),
    Smoothed_Total = slide_index_dbl(Total, Date, sum,
                                     .before = days(3), .after = days(3)),
    Smoothed_Proportion = Smoothed_count / Smoothed_Total
  ) %>%
  ungroup() %>%
  select(Date, Variant, Region, Proportion, Smoothed_Proportion)

# -- 2. colour map ----------------------------------------------------------
VariantList <- c("Alpha", "Delta", "Omicron", "BA.1", "BA.2", "BA.3", "BA.4", "BA.5",
                 "XBB", "JN.1", "XEC", "XFG", "XDV", "OtherRecombVar", "OtherVar")
ColorList   <- c("#462EB9", "#3E58CF", "#4580CA", "#549DB2", "#64BD92", "#69B091",
                 "#83BA70", "#A2BE57", "#C1BA47", "#D9AD3D", "#E69136", "#E4632E",
                 "#DC2F24", "grey70", "grey90")
Color_mapping <- setNames(ColorList, VariantList)

# -- 3. plot ----------------------------------------------------------------
p <- ggplot(plot_data, aes(x = Date, y = Smoothed_Proportion, fill = Variant)) +
  geom_area(linewidth = 0.2, colour = "black") +
  facet_wrap(~ Region, ncol = 1) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months",
               limits = c(as.Date("2020-01-01"), as.Date("2026-01-01")),
               expand = expansion(mult = 0)) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = Color_mapping) +
  labs(x = "Date", y = "Proportion", fill = "Variant") +
  theme_bw() +
  guides(fill = guide_legend(ncol = 4, byrow = FALSE)) +
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

print(p)
ggsave(p, filename = file.path(DIR_FIG_MAIN, "Figure1a_allvar_proportion_area.pdf"),
       device = cairo_pdf, width = 60, height = 125, units = "mm")
