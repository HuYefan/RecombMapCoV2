# Figure S3 — correlation matrix (VCM) of recombination variables

Pairwise correlation analysis comparing the averaged uniformed recombination
intensity with the key sequence and structural features:

* RNA fragment genomic position
* averaged GC ratio (GC%)
* averaged minimum free energy (ΔG)
* averaged paired nucleotide proportion (Proppaired)
* averaged nucleotide diversity (π)

Scatter plots and Spearman correlation coefficients are displayed for each
pairwise comparison, and the correlations are stratified by recombination
region type (**Cold**, **Warm**, **Hot**) to highlight context-dependent
dynamics. Asterisks denote significance (\* p < 0.05, \*\* p < 0.01,
\*\*\* p < 0.001).

## Script

|Script|Purpose|
|-|-|
|`scripts/01\_vcm\_correlation\_matrix.R`|`GGally::ggpairs` correlation matrix + numeric Spearman table per stratum|

## Input

|File|Description|
|-|-|
|`../data/derived/SARS2\_Results\_260811.xlsx`|**distributed** — 2,970 windows × (position, GC, ΔG, paired proportion, π, recombination intensity, TRS)|

## Notes

* Region types are the bottom / top 10% of `Recombination\_Intensity`
(`CUTOFF <- 0.10`), the same classification used for Figure 2b.
* `Genomic\_Position\_start` is used as the "RNA fragment genomic position"
variable. If you prefer the window midpoint, replace it with
`(Genomic\_Position\_start + Genomic\_Position\_end) / 2` — the manuscript's
published figure uses the window start.
* The original analysis also exported VIF, standardised regression
coefficients, an elastic-net/Lasso coefficient panel and random-forest
importance; those belong to the recombination-probability analysis rather
than to this supplementary figure and are kept in
[`../docs/SCRIPT\_MAPPING.md`](../docs/SCRIPT_MAPPING.md) for reference.

