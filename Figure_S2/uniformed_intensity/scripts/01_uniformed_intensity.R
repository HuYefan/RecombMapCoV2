# =============================================================================
#  01_uniformed_intensity.R
#  Figure S2 (panels h-j) : uniformed recombination intensity at 1-nt and
#  spike amino-acid resolution, with the ORF / spike sub-domain annotation.
# =============================================================================
#
#  Input   : data/reference/recomb_intensity_per_base.csv
#            data/reference/SARS-CoV-2_ORF.xlsx
#            data/reference/SARS-CoV-2_Spike.xlsx
#  Output  : figures/supplementary/FigureS2i_uniformed_intensity_genome.pdf
#            figures/supplementary/FigureS2j_uniformed_intensity_spike_aa.pdf
#            figures/supplementary/FigureS2j_uniformed_intensity_spike_aa.csv
#
#  The figure-1 equivalents of these two panels are produced by
#  Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R and
#  Figure_1/e_recombination_intensity/scripts/01_plot_spike_intensity.R; both
#  read the same published intensity track through load_recomb_intensity().
#
#  Run from the repository root:
#      Rscript Figure_S2/uniformed_intensity/scripts/01_uniformed_intensity.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

source("config.R")

# ---------------------------------------------------------------------------
# 1. per-base intensity (published track)
# ---------------------------------------------------------------------------
per_base <- as_tibble(load_recomb_intensity())

# ---------------------------------------------------------------------------
# 2. per-amino-acid intensity on the spike
# ---------------------------------------------------------------------------
per_aa <- per_base %>%
  filter(pos >= GENE_S["start"], pos <= GENE_S["end"]) %>%
  mutate(aa_pos = ceiling((pos - GENE_S["start"] + 1) / 3)) %>%
  group_by(aa_pos) %>%
  summarise(intensity = sum(intensity), .groups = "drop")

# ---------------------------------------------------------------------------
# 3. annotation tracks
# ---------------------------------------------------------------------------
orf_regions <- as.data.frame(read_excel(file.path(DIR_DATA_REF, "SARS-CoV-2_ORF.xlsx")))

genes_nt <- tibble(
  gene  = c("S", "RBD", "S1"),
  start = c(GENE_S["start"], GENE_RBD["start"], GENE_S1["start"]),
  end   = c(GENE_S["end"],   GENE_RBD["end"],   GENE_S1["end"])
)
aa_of <- function(nt) ceiling((nt - GENE_S["start"] + 1) / 3)
genes_aa <- tibble(gene = c("RBD", "S1"),
                   start = c(aa_of(GENE_RBD["start"]), aa_of(GENE_S1["start"])),
                   end   = c(aa_of(GENE_RBD["end"]),   aa_of(GENE_S1["end"])))

thr <- mean(per_base$intensity) + 3 * sd(per_base$intensity)

# ---------------------------------------------------------------------------
# 4. panels
# ---------------------------------------------------------------------------
p_genome <- ggplot() +
  geom_rect(data = genes_nt,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_line(data = per_base, aes(x = pos, y = intensity),
            color = "#17406D", linewidth = 0.2) +
  geom_hline(yintercept = thr, color = "darkred", linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 29903), breaks = c(0, seq(0, 5) * 5000, 29903)) +
  labs(x = "Genomic Position (nt)", y = "Uniformed\n Event Count") +
  theme_ms()

p_spike <- ggplot() +
  geom_rect(data = genes_aa,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_line(data = per_aa, aes(x = aa_pos, y = intensity),
            color = "#17406D", linewidth = 0.4) +
  geom_hline(yintercept = thr, color = "darkred", linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(limits = c(1, 1274),
                     breaks = c(1, genes_aa$start[1], genes_aa$end[1],
                                genes_aa$end[2], 1274)) +
  labs(x = "Spike amino-acid position", y = "Uniformed\n Event Count\n (Per AA)") +
  theme_ms()

print(p_genome)
print(p_spike)

# the published uniformed-intensity panels are 160 x 30 mm strips
ggsave(p_genome,
       filename = file.path(DIR_FIG_SUPP, "FigureS2i_uniformed_intensity_genome.pdf"),
       device = cairo_pdf, width = 160, height = 30, units = "mm")
ggsave(p_spike,
       filename = file.path(DIR_FIG_SUPP, "FigureS2j_uniformed_intensity_spike_aa.pdf"),
       device = cairo_pdf, width = 160, height = 30, units = "mm")

aa_out <- file.path(DIR_FIG_SUPP, "FigureS2j_uniformed_intensity_spike_aa.csv")
write_csv(per_aa, aa_out)
message("wrote ", aa_out)
