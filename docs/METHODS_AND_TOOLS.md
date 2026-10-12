# Methods and upstream tools

The scripts in this repository cover the **downstream** analysis and figure
generation. The upstream steps that produce the primary inputs are listed here
with the version/parameter actually used, so the pipeline can be re-run end to
end.

## Upstream pipeline

|Step|Tool|Purpose|Output consumed by|
|-|-|-|-|
|Sequence curation|Metadata export|collection date, location, Pango lineage, AA substitutions|Figure 1a, Figure S1|
|Tree / ancestral reconstruction|[UShER](https://github.com/yatisht/usher) + [matUtils](https://github.com/yatisht/usher) on the UCSC `public-2025-12-31` dataset (8,433,953 variants)|mutation-annotated tree with recombination nodes|Figure 1c–e, Figure S2|
|Recombination calling|[RIPPLES](https://github.com/PoonLab/ripples)|per-node recombination breakpoint intervals + `recomb\_parsimony` support|Figure 1b–e, Figure S2|
|Descendant extraction|[UShER](https://github.com/yatisht/usher) / [matUtils](https://github.com/yatisht/usher) tree traversal|node → descendant strain list with collection dates|Figure 1c–e, Figure S2|
|Consensus genome|42,911-sequence consensus|reference for the per-site features|Figure 2, Figure S3|
|RNA secondary structure|[LinearFold](https://github.com/LinearFold/LinearFold)|per-sequence GC, ΔG and dot-bracket structure|Figure 2, Figure S3|
|Nucleotide diversity|per-site π and entropy across the consensus set|π track|Figure 2, Figure S3|
|Global phylogeny|[Nextstrain ncov](https://github.com/nextstrain/ncov) "open / global / all-time" build|subsampled tree JSON + tip metadata|Figure S1|
|Recombination pathway diagrams|[MAFFT](https://mafft.cbrc.jp/alignment/software/) + [snipit](https://github.com/aineniamh/snipit)|per-lineage recombination-mode alignment plots|Figure 1b|
|Spike structure|PDB [6VXX](https://www.rcsb.org/structure/6VXX)|WT spike ectodomain trimer|Figure 1f–g|

## Thresholds and parameters

|Parameter|Value|Where|
|-|-|-|
|Recombination confidence filter|`recomb\_parsimony >= 7`|`config.R` → `RECOMB\_PARSIMONY\_MIN`|
|Breakpoint interval weighting|each interval of length *L* contributes 1/*L* per base|`config.R` → `uniformed\_intensity()`|
|Intensity track used|`data/reference/recomb\_intensity\_per\_base.csv`, read via `config.R` → `load\_recomb\_intensity()`|`config.R`|
|Sliding window (features)|100 nt window, 10 nt step|`Figure\_2/scripts/02\_feature\_windows.R`|
|Cold / hot regions|bottom / top 10% of recombination intensity|`Figure\_2/scripts/03\_plot\_genome\_features.R`, `Figure\_S3/scripts/01\_vcm\_correlation\_matrix.R`|
|Permutation test|`perm::permTS`, `exact.mc`, N = 5,000, one-sided "greater"|`Figure\_S2/permutation\_hotspots/scripts/01\_permutation\_test.R`|
|Multiple testing|Benjamini-Hochberg (FDR < 0.05)|permutation test, correlation table|
|Correlation method|Spearman rank, `exact = FALSE`|`Figure\_S3/scripts/01\_vcm\_correlation\_matrix.R`|
|Cut-off search (supporting)|5–50% in 5% steps, Wilcoxon hot vs. cold|documented in `docs/SCRIPT\_MAPPING.md`|

## Genome coordinates (Wuhan-Hu-1, NC\_045512.2)

All coordinates live in `config.R` / `config.py`:

|Feature|Start|End|
|-|-|-|
|Genome|1|29,903|
|Spike (S) CDS|21,563|25,384|
|S1 subunit|21,563|23,617|
|RBD|22,517|23,185|
|RBM (hottest breakpoint window)|22,895|22,927|

## Data layout

```
data/raw/         primary inputs the scripts read (recombination events,
                  aggregated variant proportions, node dates)
data/reference/   annotated inputs (uniformed intensity track, ORF coordinates,
                  variant mappings, 6VXX, Nextstrain tree + metadata)
data/derived/     per-site / windowed feature tables (SARS2\_Results.xlsx,
                  regenerated intermediates)
```

Only files that directly serve Figures 1–2 and Supplementary Figures S2–S3 are
kept; see [`../data/README.md`](../data/README.md) for the full list and for the
large primary inputs that are not distributed.

