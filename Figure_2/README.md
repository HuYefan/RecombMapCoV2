# Figure 2 — genomic features of recombination hotspots

**Panel a.** Genome-wide distribution of the averaged sequence and structural
features mapped alongside the averaged uniformed recombination intensity
(bottom track): averaged GC ratio (GC%), averaged minimum free energy (ΔG),
averaged paired nucleotide proportion (Proppaired) and averaged nucleotide
diversity (π).

**Panel b.** Cold regions (bottom 10% of recombination intensity) vs. hot
regions (top 10%), compared with two-sided Mann-Whitney tests.

## Scripts

|Script|Purpose|
|-|-|
|`scripts/01\_feature\_summary.R`|collapse the per-sequence RNA structure predictions into per-position means|
|`scripts/02\_feature\_windows.R`|100-nt / 10-nt step sliding-window table joining all features|
|`scripts/03\_plot\_genome\_features.R`|panel a (feature tracks) + panel b (boxplots) + statistics table|

## Inputs

|File|Description|
|-|-|
|`../data/derived/SARS2\_Results\_260811.xlsx`|**distributed** — windowed feature table (2,970 windows × 10 columns)|
|`../data/reference/recomb\_intensity\_per\_base.csv`|1-nt uniformed intensity (from Figure 1c–d)|
|`../data/raw/SARS2\_Structures\_Detailed.csv`|\~300 MB per-sequence folding table — not distributed|
|`../data/derived/SARS2\_Genome\_Diversity\_Entropy.csv`|per-position π and entropy — not distributed|

Because `SARS2\_Results\_260811.xlsx` is distributed, `03\_plot\_genome\_features.R`
runs out of the box. Steps 01 and 02 are only needed to rebuild that table from
the raw folding / diversity outputs.

