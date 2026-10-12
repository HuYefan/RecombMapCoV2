# =============================================================================
#  01_plot_spike_intensity.R
#  Figure 1e : high-resolution mapping of the uniformed recombination
#  intensity across the spike gene, at amino-acid resolution.
# =============================================================================
#
#  Input   : data/reference/recomb_intensity_per_base.csv
#            (published per-base uniformed intensity track, read through
#             load_recomb_intensity() in config.R)
#  Output  : figures/main/Figure1e_spike_intensity_aa.pdf
#            figures/main/Figure1e_spike_intensity_aa.csv
#            (per-residue intensity, reused by the structural mapping of
#             Figure 1f-g and by Figure S2j)
#
#  Method
#  ------
#  The 1-nt uniformed intensity is restricted to the spike CDS
#  (21563-25384, NC_045512.2) and summed over the three nucleotides of each
#  codon to give a per-amino-acid intensity (S1-S1274).  The dashed line marks
#  the mean + 3 SD of the genome-wide per-base intensity.
#
#  Run from the repository root:
#      Rscript Figure_1/e_recombination_intensity/scripts/01_plot_spike_intensity.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("config.R")

recomb_results_nt <- as_tibble(load_recomb_intensity())

# -- convert nt intensity to amino-acid intensity ----------------------------
recomb_results_aa <- recomb_results_nt %>%
  filter(pos >= GENE_S["start"], pos <= GENE_S["end"]) %>%
  mutate(aa_pos = ceiling((pos - GENE_S["start"] + 1) / 3)) %>%
  group_by(aa_pos) %>%
  summarise(intensity = sum(intensity), .groups = "drop")

# -- annotation bands in amino-acid coordinates ------------------------------
aa_of <- function(nt) ceiling((nt - GENE_S["start"] + 1) / 3)
genes_aa <- tibble(
  gene  = c("RBD", "S1"),
  start = c(aa_of(GENE_RBD["start"]), aa_of(GENE_S1["start"])),
  end   = c(aa_of(GENE_RBD["end"]),   aa_of(GENE_S1["end"]))
)

# threshold defined on the genome-wide per-base intensity
thr <- mean(recomb_results_nt$intensity) + 3 * sd(recomb_results_nt$intensity)

p <- ggplot() +
  geom_rect(data = genes_aa,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.2, fill = "#0A7ECC") +
  geom_line(data = recomb_results_aa, aes(x = aa_pos, y = intensity),
            color = "#17406D", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(limits = c(1, 1274),
                     breaks = c(1, genes_aa$start[1], genes_aa$end[1],
                                genes_aa$end[2], 1274)) +
  geom_hline(yintercept = thr, color = "darkred", linewidth = 0.5, linetype = "dashed") +
  labs(x = "Spike amino-acid position",
       y = "Uniformed\n Event Count\n (Per AA)") +
  theme_ms()

print(p)
ggsave(p, filename = file.path(DIR_FIG_MAIN, "Figure1e_spike_intensity_aa.pdf"),
       device = cairo_pdf, width = 160, height = 45, units = "mm")

# per-residue table next to the figure (consumed by Figure 1f-g and Figure S2j)
aa_out <- file.path(DIR_FIG_MAIN, "Figure1e_spike_intensity_aa.csv")
write_csv(recomb_results_aa, aa_out)
message("wrote ", aa_out)
