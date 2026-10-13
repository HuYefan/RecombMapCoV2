# Data

Only inputs that directly produce **Figures 1–2 and Supplementary Figures S1–S3**
are kept here.

```
data/
├── raw/            primary inputs the scripts read
├── reference/      annotated inputs (coordinates, structures, variant lists)
└── derived/        per-site / windowed feature tables
```

## `raw/`

|File|Size|Used by|
|-|-|-|
|`hcov-19_proportion.csv`|1.8 MB|Figure 1a — variant proportions per collection date and location (recombinant groups)|
|`hcov-19_VOIproportion.csv`|3.0 MB|Figure 1a (all-variant panel) — variant proportions per date and location|
|`recombination.tsv`|9.0 MB|Figures 1c–e, 2a, S2 — RIPPLES recombination events, 1,273 nodes with `recomb_parsimony >= 7`|
|`descendants_node_dates.tsv`|47 KB|Figure S2 — earliest descendant collection date per recombination node|

## `reference/`

|File|Size|Used by|
|-|-|-|
|`recomb_intensity_per_base.csv`|541 KB|Figures 1c–e, 2a, S2i–j — per-base uniformed recombination intensity (`pos`, `intensity`; positions 0–29 902)|
|`SARS-CoV-2_ORF.xlsx` / `.csv`|15 KB|Figures 1c–e, S2 — ORF coordinates (NC_045512.2)|
|`SARS-CoV-2_Spike.xlsx`|15 KB|Figure S2 — spike sub-domain coordinates|
|`RecombVar_ofInterest.xlsx`|16 KB|Figure 1a — recombinant group → sublineage prefixes|
|`Variants_of_Interest.xlsx`|16 KB|Figure 1a (all-variant panel) — variant → sublineage prefixes|
|`Recombinant_Variants.xlsx`|206 KB|Figure S2 — designated recombinants and their breakpoint coordinates|
|`6vxx.pdb`|4.7 MB|Figure 1f–g — WT spike ectodomain trimer (PDB 6VXX)|
|`ncov_tree.json`|4.2 MB|Figure S1 — Nextstrain v2 tree, 5,714 nodes / 2,938 tips|
|`ncov_metadata.tsv`|658 KB|Figure S1 — matching tip metadata, complete 2,938-tip set|


## `derived/`

|File|Size|Used by|
|-|-|-|
|`SARS2_Results.xlsx`|225 KB|Figures 2 and S3 — windowed feature table, 2,970 windows × (position start/end, Mean_GC, Mean_Energy, Paired_Proportion, Unpaired_Proportion, Pi_Site, Recombination_Intensity, Sample_Count, TRS)|
|`SARS2_summary_results.csv`|—|*(regenerated)* Figure 2 step 1 output — per-position Mean_GC / Mean_Energy / Dot_Proportion|
|`SARS2_Genome_Diversity_Entropy.csv`|—|*(not distributed)* per-position π and entropy|
|`SARS2_Feature_Windows.csv`|—|*(regenerated)* Figure 2 step 2 output — 100-nt / 10-nt stepped feature windows|
|`recombination_events_summary.tsv`|—|*(regenerated)* one row per event, from `annotate_recombination_dates.py`|

## Files that must be obtained separately

|File|Approx. size|Goes to|Source|
|-|-|-|-|
|`SARS2_Structures_Detailed.csv`|300 MB|`data/raw/`|per-sequence RNA secondary-structure predictions (LinearFold) for 42,911 consensus genomes|
|`SARS2_Genome_Diversity_Entropy.csv`|455 KB|`data/derived/`|per-position nucleotide diversity / entropy across the same consensus set|
|`descendants.tsv`|669 MB|`data/raw/`|UShER/RIPPLES node → descendant sequence list, needed only to regenerate `descendants_node_dates.tsv`|
|`NC_045512.2.fasta` / `.gb`|30 / 80 KB|`data/reference/`|[NCBI NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) reference genome|
|`<LINEAGE>.fasta`|—|`Figure_1/b_recombination_pathways/fasta/`|GISAID sequences for the Figure 1b pathway diagrams — accessions are listed in the script and in that folder's README|

## Recreating `recombination.tsv`

`recombination.tsv` is RIPPLES output: for every candidate recombination node it
reports the two breakpoint intervals, the donor and acceptor node, and the
`recomb_parsimony` support score. All figures apply `recomb_parsimony >= 7`
(`RECOMB_PARSIMONY_MIN` in `config.R`). The caller is not part of this
repository; see [`../docs/METHODS\_AND\_TOOLS.md`](../docs/METHODS_AND_TOOLS.md).

