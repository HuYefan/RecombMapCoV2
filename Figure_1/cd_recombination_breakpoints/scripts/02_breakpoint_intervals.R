# =============================================================================
#  02_breakpoint_intervals.R
#  Figure 1c-d and Figure S2 : breakpoint interval map, midpoint histogram and
#  the per-base uniformed recombination intensity track.
# =============================================================================
#
#  Input   : data/raw/recombination.tsv                      (RIPPLES events)
#  Output  : figures/main/Figure1c_breakpoint_intervals.pdf
#            figures/main/Figure1d_breakpoint_midpoint_histogram.pdf
#            figures/main/Figure1e_uniformed_intensity_genome.pdf
#            data/reference/recomb_intensity_per_base.csv   (shared with Fig. 2 / S2)
#
#  Method
#  ------
#  Each recombination event contributes uniformly to every base of each of its
#  two breakpoint intervals, i.e. an interval [start, end] of length L gives a
#  weight of 1/L to each of its L bases.  Summing over all events gives the
#  "uniformed recombination intensity" at 1-nt resolution.
#
#  Run from the repository root:
#      Rscript Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("config.R")

# -- 1. read and filter the recombination calls ------------------------------
df <- read_tsv(file.path(DIR_DATA_RAW, "recombination.tsv"), show_col_types = FALSE)

df_filtered <- df %>%
  filter(recomb_parsimony >= RECOMB_PARSIMONY_MIN) %>%
  mutate(
    `breakpoint-1_interval` = str_replace_all(`breakpoint-1_interval`, "GENOME_SIZE", "29903"),
    `breakpoint-2_interval` = str_replace_all(`breakpoint-2_interval`, "GENOME_SIZE", "29903")
  )

# -- 2. parse "(start,end)" intervals ----------------------------------------
extract_coords <- function(data, col_name, prefix) {
  data %>%
    tidyr::extract(
      col = !!sym(col_name),
      into = c(paste0(prefix, "_start"), paste0(prefix, "_end")),
      regex = "\\((\\d+),(\\d+)\\)",
      convert = TRUE,
      remove = FALSE
    )
}

plot_data <- df_filtered %>%
  extract_coords("breakpoint-1_interval", "bp1") %>%
  extract_coords("breakpoint-2_interval", "bp2") %>%
  mutate(row_id = row_number()) %>%
  filter(bp1_start != 0,
         bp1_end - bp1_start > 1,
         bp2_end != 29903,
         bp2_end - bp2_start > 1)

message("breakpoint intervals retained: ", nrow(plot_data))

# -- 3. panel d: midpoint histogram data -------------------------------------
bp1_mid <- plot_data %>%
  filter(!is.na(bp1_start)) %>%
  mutate(mid = (bp1_start + bp1_end) / 2) %>%
  select(mid)

bp2_mid <- plot_data %>%
  filter(!is.na(bp2_start)) %>%
  mutate(mid = (bp2_start + bp2_end) / 2) %>%
  select(mid)

density_data <- bind_rows(bp1_mid, bp2_mid)

# -- 4. per-base uniformed intensity ----------------------------------------
parse_and_expand <- function(df) {
  df %>%
    mutate(
      s1 = as.numeric(str_extract(`breakpoint-1_interval`, "(?<=\\()[0-9]+")),
      e1 = as.numeric(str_extract(`breakpoint-1_interval`, "[0-9]+(?=\\))")),
      s2 = as.numeric(str_extract(`breakpoint-2_interval`, "(?<=\\()[0-9]+")),
      e2 = as.numeric(str_extract(`breakpoint-2_interval`, "[0-9]+(?=\\))"))
    ) %>%
    pivot_longer(cols = starts_with("s"), names_to = "type", values_to = "start") %>%
    mutate(end = if_else(type == "s1", e1, e2)) %>%
    select(start, end) %>%
    filter(!is.na(start))
}

intervals <- parse_and_expand(df_filtered) %>%
  mutate(
    length = end - start,
    weight = if_else(length > 0, 1 / length, 1)
  )

# -- 4b. per-base uniformed intensity ----------------------------------------
# Read through the shared helper so Figure 1e, Figure 2a and Figure S2 use
# exactly the same values (data/reference/recomb_intensity_per_base.csv).
recomb_results <- as_tibble(load_recomb_intensity())

# -- 5. annotation rectangles (spike / RBD / S1) -----------------------------
genes <- tibble(
  gene  = c("S", "RBD", "S1"),
  start = c(GENE_S["start"], GENE_RBD["start"], GENE_S1["start"]),
  end   = c(GENE_S["end"],   GENE_RBD["end"],   GENE_S1["end"])
)
x_breaks <- c(0, seq(0, 5) * 5000, 29903)

# -- panel c : breakpoint interval swimlane ----------------------------------
plot_swimlane <- plot_data %>%
  arrange(desc(bp2_end), bp2_start) %>%
  mutate(y_position = row_number())

p_c <- ggplot(plot_swimlane) +
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_rect(aes(xmin = bp1_start, xmax = bp1_end,
                ymin = y_position - 0.4, ymax = y_position + 0.4),
            fill = "#0A7ECC") +
  geom_rect(data = plot_swimlane %>% filter(!is.na(bp2_start)),
            aes(xmin = bp2_start, xmax = bp2_end,
                ymin = y_position - 0.4, ymax = y_position + 0.4),
            fill = "#0297D6") +
  scale_x_continuous(limits = c(0, 29903), breaks = x_breaks) +
  labs(x = "Genomic Position (nt)", y = "Recombination\n Events") +
  theme_ms() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

# -- panel d : midpoint histogram --------------------------------------------
p_d <- ggplot(density_data, aes(x = mid)) +
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_histogram(binwidth = 1, fill = "#4DBBD5", color = "#17406D", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 29903), breaks = x_breaks) +
  labs(x = "Genomic Position (nt)", y = "Recombination\n Midpoint\n Count") +
  theme_ms()

# -- panel e : uniformed intensity across the genome -------------------------
p_e <- ggplot() +
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_line(data = recomb_results, aes(x = pos, y = intensity),
            color = "#17406D", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 29903), breaks = x_breaks) +
  labs(x = "Genomic Position (nt)", y = "Uniformed\n Recombination\n Event Count") +
  theme_ms()

# -- save --------------------------------------------------------------------
gsave_pdf <- function(plot, file, width, height) {
  ggsave(filename = file.path(DIR_FIG_MAIN, file), plot = plot,
         device = cairo_pdf, width = width, height = height, units = "mm")
  message("wrote ", file.path(DIR_FIG_MAIN, file))
}

gsave_pdf(p_c, "Figure1c_breakpoint_intervals.pdf", 160, 90)
gsave_pdf(p_d, "Figure1d_breakpoint_midpoint_histogram.pdf", 160, 50)
# the published uniformed-intensity panel is a 160 x 30 mm strip
gsave_pdf(p_e, "Figure1e_uniformed_intensity_genome.pdf", 160, 30)

# combined c/d/e strip used when assembling the manuscript figure
final_plot <- p_c / p_d / p_e +
  plot_layout(heights = c(3, 1, 1)) &
  theme(axis.title = element_text(size = 7), title = element_text(size = 7))
ggsave(final_plot, filename = file.path(DIR_FIG_MAIN, "Figure1cde_breakpoint_strip.pdf"),
       device = cairo_pdf, width = 160, height = 120, units = "mm")

print(final_plot)
