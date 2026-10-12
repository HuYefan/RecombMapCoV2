# Figure 1c–d — recombination breakpoint intervals and midpoint histogram

**Panel c.** Genome-wide mapping of the 2,991 recombination breakpoint intervals
derived from 8,433,953 UCSC variants. Each interval is drawn as a horizontal
segment; the spike, S1 and RBD regions are shaded. Confidence filter:
`recomb_parsimony >= 7`.

**Panel d.** Histogram of breakpoint midpoints (bin = 1 nt). Recombination
events can only be localised to the interval between two recombination-
informative SNPs, so the midpoint of each prediction interval is used.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/01_annotate_recombination_dates.py` | attach the earliest descendant collection date to every recombination event (reads the full `descendants.tsv`) |
| `scripts/02_breakpoint_intervals.R` | swimlane plot of all breakpoint intervals (panel c), midpoint histogram + per-base uniformed intensity (panel d) |
| `scripts/03_plot_events_over_time.R` | monthly counts / cumulative number of new recombination nodes (supporting analysis) |

## Inputs

| File | Description |
|---|---|
| `../../data/raw/recombination.tsv` | RIPPLES recombination events (one row per candidate breakpoint) |
| `../../data/raw/descendants_node_dates.tsv` | compact per-node earliest-date table (1273 nodes) used by the R code |
| `../../data/reference/SARS-CoV-2_ORF.xlsx` | ORF coordinates for the per-region views |

The upstream `descendants.tsv` (669 MB) is not distributed. Put it in
`data/raw/` and run `scripts/01_annotate_recombination_dates.py` to regenerate
`recombination_annotated_with_dates.tsv`, `recombination_events_summary.tsv`
and the compact `descendants_node_dates.tsv` table.

## Outputs

* `../../figures/main/Figure1c_breakpoint_intervals.pdf`
* `../../figures/main/Figure1d_breakpoint_midpoint_histogram.pdf`
* `../../figures/main/Figure1e_uniformed_intensity_genome.pdf`
* `../../figures/main/Figure1cde_breakpoint_strip.pdf` (c/d/e assembled, matches the published `Recomb_Distribution_260207.pdf`)

## Notes

The **uniformed recombination intensity** — each breakpoint interval of length
`L` contributing weight `1/L` to every base it spans — is a shared input here,
in Figure 1e, in the bottom track of Figure 2a and in Figure S2i–j. It is read
from `../../data/reference/recomb_intensity_per_base.csv` via
`config.R::load_recomb_intensity()`, so every panel uses identical values; the
re-derivation is `config.R::uniformed_intensity()`, covered by
`../../tests/test_uniformed_intensity.R`.
