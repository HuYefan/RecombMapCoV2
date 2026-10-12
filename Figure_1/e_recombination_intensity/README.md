# Figure 1e — uniformed recombination intensity along the genome and the spike

High-resolution mapping of the uniformed recombination intensity, highlighting
the hotspot regions inside the spike gene. This is the same quantity shown in
Figure 1c–d (bottom track), Figure 2a (bottom track) and Figure S2i–j, at
1-nt and spike amino-acid resolution.

## Script

| Script | Purpose |
|---|---|
| `scripts/01_plot_spike_intensity.R` | restrict the 1-nt intensity to the spike CDS, sum per codon → per-amino-acid intensity plot |

## Input

| File | Description |
|---|---|
| `../../data/reference/recomb_intensity_per_base.csv` | 1-nt uniformed intensity — the track shared with Figures 1c–d, 2a and S2i–j, read via `config.R::load_recomb_intensity()` |

## Outputs

* `../../figures/main/Figure1e_spike_intensity_aa.pdf`
* `../../figures/main/Figure1e_spike_intensity_aa.csv` (per-residue intensity, reused by the structural mapping of Figure 1f–g and Figure S2j)

## Notes

* Spike CDS coordinates (NC_045512.2): 21563–25384; S1: 21563–23617; RBD:
  22517–23185; RBM hotspot: 22895–22927 (S:445–455).
* The dashed threshold is mean + 3 SD of the **genome-wide** per-base intensity.
