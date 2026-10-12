# Figure S2 — detailed nucleotide breakpoints and genomic recombination intensity

**Panels a–g.** Recombination pathways and per-lineage breakpoint intervals for
XBB, XDE, XDV, XEC and XFG, with the corresponding parental crosses.

**Panels h–j.** Detailed nucleotide breakpoints, the per-base uniformed
recombination intensity across the genome and along the spike gene, and the
Monte-Carlo permutation test for regions with significantly elevated
recombination intensity.

Panels a–g reuse the plotters of Figure 1b/1c–e on the per-lineage inputs; the
permutation test lives here.

## Scripts

|Script|Purpose|
|-|-|
|`permutation\_hotspots/scripts/01\_permutation\_test.R`|`perm::permTS` Monte-Carlo test (N = 5,000) per genomic region, BH-corrected|
|`breakpoint\_intervals/scripts/01\_breakpoint\_intervals.R`|per-lineage swimlane of breakpoint intervals (panels a–g) — wrapper, see notes|
|`uniformed\_intensity/scripts/01\_uniformed\_intensity.R`|per-base and per-amino-acid uniformed intensity tracks (panels h–j)|

## Shared code (not duplicated)

The heavy lifting for panels a–g and h–j already exists once in this repository;
run those scripts and reuse their outputs:

|What|Where|
|-|-|
|breakpoint interval map, midpoint histogram, per-base intensity|[`../Figure\_1/cd\_recombination\_breakpoints/scripts/02\_breakpoint\_intervals.R`](../Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R)|
|spike amino-acid intensity|[`../Figure\_1/e\_recombination\_intensity/scripts/01\_plot\_spike\_intensity.R`](../Figure_1/e_recombination_intensity/scripts/01_plot_spike_intensity.R)|
|per-lineage recombination-pathway diagrams|[`../Figure\_1/b\_recombination\_pathways/recombination\_pathways.py`](../Figure_1/b_recombination_pathways/recombination_pathways.py)|
|event dating (drives panel a–g ordering)|[`../Figure\_1/cd\_recombination\_breakpoints/scripts/01\_annotate\_recombination\_dates.py`](../Figure_1/cd_recombination_breakpoints/scripts/01_annotate_recombination_dates.py)|

`breakpoint\_intervals/scripts/01\_breakpoint\_intervals.R` and
`uniformed\_intensity/scripts/01\_uniformed\_intensity.R` are thin wrappers that
call those same functions with the Figure S2 output names, so the supplementary
figure can be regenerated on its own.

## Inputs

|File|Description|
|-|-|
|`../data/raw/recombination.tsv`|**distributed** — RIPPLES recombination events|
|`../data/raw/descendants\_node\_dates.tsv`|**distributed** — per-node earliest collection date|
|`../data/reference/SARS-CoV-2\_ORF.xlsx`|**distributed** — ORF coordinates|
|`../data/reference/SARS-CoV-2\_Spike.xlsx`|**distributed** — spike sub-domain coordinates|



