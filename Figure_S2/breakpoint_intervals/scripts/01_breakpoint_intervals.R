# =============================================================================
#  01_breakpoint_intervals.R
#  Figure S2 (panels a-g) : per-lineage breakpoint interval map, reusing the
#  Figure 1c breakpoint parser on the same underlying event table.
# =============================================================================
#
#  The parsing, weighting and plotting logic is identical to
#  Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R.
#  This wrapper re-runs it with the Figure S2 output names and a per-lineage
#  breakpoint panel built from the designation table.
#
#  Input   : data/raw/recombination.tsv
#            data/reference/Recombinant_Variants.xlsx  (per-lineage breakpoints)
#            data/reference/SARS-CoV-2_ORF.xlsx
#  Output  : figures/supplementary/FigureS2_breakpoint_intervals.pdf
#            figures/supplementary/FigureS2_lineage_breakpoints.pdf
#
#  Run from the repository root:
#      Rscript Figure_S2/breakpoint_intervals/scripts/01_breakpoint_intervals.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})

source("config.R")

# ---------------------------------------------------------------------------
# 1. events -> breakpoint intervals (identical logic to Figure 1c-d)
# ---------------------------------------------------------------------------
df <- read_tsv(file.path(DIR_DATA_RAW, "recombination.tsv"), show_col_types = FALSE)

df_filtered <- df %>%
  filter(recomb_parsimony >= RECOMB_PARSIMONY_MIN) %>%
  mutate(
    `breakpoint-1_interval` = str_replace_all(`breakpoint-1_interval`, "GENOME_SIZE", "29903"),
    `breakpoint-2_interval` = str_replace_all(`breakpoint-2_interval`, "GENOME_SIZE", "29903")
  )

extract_coords <- function(data, col_name, prefix) {
  data %>%
    tidyr::extract(col = !!sym(col_name),
                   into = c(paste0(prefix, "_start"), paste0(prefix, "_end")),
                   regex = "\\((\\d+),(\\d+)\\)", convert = TRUE, remove = FALSE)
}

plot_data <- df_filtered %>%
  extract_coords("breakpoint-1_interval", "bp1") %>%
  extract_coords("breakpoint-2_interval", "bp2") %>%
  mutate(row_id = row_number()) %>%
  filter(bp1_start != 0, bp1_end - bp1_start > 1,
         bp2_end != 29903, bp2_end - bp2_start > 1)

genes <- tibble(
  gene  = c("S", "RBD", "S1"),
  start = c(GENE_S["start"], GENE_RBD["start"], GENE_S1["start"]),
  end   = c(GENE_S["end"],   GENE_RBD["end"],   GENE_S1["end"])
)
x_breaks <- c(0, seq(0, 5) * 5000, 29903)

p_events <- plot_data %>%
  arrange(desc(bp2_end), bp2_start) %>%
  mutate(y_position = row_number()) %>%
  ggplot() +
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_rect(aes(xmin = bp1_start, xmax = bp1_end,
                ymin = y_position - 0.4, ymax = y_position + 0.4),
            fill = "#0A7ECC") +
  scale_x_continuous(limits = c(0, 29903), breaks = x_breaks) +
  labs(x = "Genomic Position (nt)", y = "Recombination\n Events") +
  theme_ms() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

print(p_events)
ggsave(p_events,
       filename = file.path(DIR_FIG_SUPP, "FigureS2_breakpoint_intervals.pdf"),
       device = cairo_pdf, width = 160, height = 90, units = "mm")

# ---------------------------------------------------------------------------
# 2. designated-recombinant breakpoints per lineage (the 173 designated strains)
# ---------------------------------------------------------------------------
designation_file <- file.path(DIR_DATA_REF, "Recombinant_Variants.xlsx")
if (file.exists(designation_file)) {
  recombination <- as.data.frame(read_excel(designation_file, col_names = TRUE))

  recombination_final <- recombination %>%
    separate_rows(Breakpoint, sep = ",") %>%
    filter(Breakpoint != "") %>%
    separate(Breakpoint, into = c("BP_start", "BP_end"), sep = "-", convert = TRUE) %>%
    mutate(across(c(BP_start, BP_end), as.numeric))

  recombination_final$Lineage <- factor(recombination_final$Lineage,
                                        levels = rev(unique(recombination$Lineage)))

  p_lineage <- ggplot(recombination_final,
                      aes(x = BP_start, xend = BP_end, y = Lineage, yend = Lineage)) +
    geom_segment(color = "#4DBBD5", linewidth = 1) +
    geom_rect(data = genes,
              aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
    scale_x_continuous(limits = c(1, 29903),
                       breaks = c(1, 5000, 10000, 15000, 20000, 25000, 29903)) +
    labs(x = "Genomic Position (nt)", y = "Lineage") +
    theme_ms() +
    theme(panel.grid.major.x = element_line(color = "grey80"),
          panel.grid.major.y = element_blank())

  print(p_lineage)
  ggsave(p_lineage,
         filename = file.path(DIR_FIG_SUPP, "FigureS2_lineage_breakpoints.pdf"),
         device = cairo_pdf, width = 160, height = 120, units = "mm")
} else {
  message("[skip] ", designation_file, " not found; per-lineage panel not drawn")
}
