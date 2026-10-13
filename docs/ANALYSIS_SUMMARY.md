# Analysis summary — algorithms and data generated

This document summarises **what the code computes** (the algorithms) and **what
data it produces** (the generated datasets), figure by figure. It is written to
be usable directly when filling in a Data Availability statement.

Repository: **https://github.com/HuYefan/RecombMapCoV2** (code `RecombMapCoV2`)


## 1\. Algorithms

### 1.1 Variant prevalence (Figure 1a)

* **Lineage grouping.** Each Pango lineage is reduced to its prefix (text before
the first `.`), `(consensus call)` suffixes are stripped, and the prefix is
mapped to a recombinant group through the `Sublineage` field of
`RecombVar_ofInterest.xlsx` / `Variants_of_Interest.xlsx` (a value such as
`XDV|NB` expands to several prefixes). Unmatched lineages become `Others`.
* **Proportions.** Counts are aggregated per collection date and per location
(country), and a global `World` aggregate is appended
(`Proportion = count / Total`).
* **Smoothing.** A **7-day centred rolling sum** is applied to both numerator and
denominator before the ratio is taken
(`Smoothed_Proportion = Smoothed_count / Smoothed_Total`), which avoids gaps
and unstable ratios on low-count days.
* **Produced by.** `Figure_1/a_variant_proportion/scripts/01`–`05`.

### 1.2 Recombination-pathway delineation (Figure 1b)

* For each recombinant lineage the lineage-defining sequences and the two
parental references are MAFFT-aligned (`--auto`) and drawn in snipit
**recombination mode** (`--recombi-mode --recombi-references`), which assigns
every alignment column to the parental lineage it most resembles and infers the
switch points.
* **Produced by.** `Figure_1/b_recombination_pathways/recombination_pathways.py`.

### 1.3 Uniformed recombination intensity (Figures 1c–1e, 2a, S2)

* **Event filter.** Only recombination nodes with `recomb_parsimony >= 7` are
used (1,273 nodes in RIPPLES output; 55,498 breakpoint intervals).
* **Interval weighting.** A recombination event localises each breakpoint to an
interval between two recombination-informative SNPs, so the event is spread
**uniformly** over that interval: an interval of length *L* contributes weight
**1/L** to each of its *L* bases. The per-base **uniformed recombination
intensity** is the sum over all intervals covering a base.
* **Implementation.** A difference array with one cumulative sum, with the
weights aggregated per index (`rowsum`) — an O(n) equivalent of the per-base
loop. In R a naive `delta[idx] <- delta[idx] + weight` is not an accumulating
update when the subscript repeats, so the weights must be summed per index
first. The track the submitted figures were drawn from ships as
`data/reference/recomb_intensity_per_base.csv` and is what every panel reads
through `config.R::load_recomb_intensity()`; the re-derivation is kept in
`config.R::uniformed_intensity()` and covered by
`tests/test_uniformed_intensity.R`.
* **Peak regions.** nt ≈ 284–302 and the RBM window nt 22,895–22,927 (S:445–455).
* **Amino-acid resolution.** The 1-nt intensity restricted to the spike CDS
(21,563–25,384) is summed over each codon, giving a per-residue intensity for
S1–S1274.
* **Threshold.** mean + 3 SD of the genome-wide per-base intensity (dashed line
in the figures).
* **Peak regions (both tracks).** nt ≈ 284–302 and the RBM window
nt 22,895–22,927 (S:445–455).
* **Regression test.** `tests/test_uniformed_intensity.R` (8 checks, including a
comparison against the straightforward per-base reference loop on the real
event table).
* **Produced by.** `Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R`,
`Figure_1/e_recombination_intensity/scripts/01_plot_spike_intensity.R`,
`Figure_S2/uniformed_intensity/scripts/01_uniformed_intensity.R`.

### 1.4 Per-base midpoint histogram (Figure 1d)

* Because an event can only be localised to an interval, the **midpoint** of each
breakpoint interval is used as the point estimate of the breakpoint position;
the histogram uses a **bin width of 1 nt**.
* **Produced by.** `Figure_1/cd_recombination_breakpoints/scripts/02_breakpoint_intervals.R`.

### 1.5 Event dating (Figure S2, supporting panels)

* Each recombination node is dated by the **earliest collection date among its
descendant sequences** (an upper bound on the true recombination time). Only
complete `YYYY-MM-DD` dates are used for the primary estimate; nodes without
one fall back to partial dates (`YYYY-MM-01`, `YYYY-01-01`) and are flagged
`partial_only`.
* **Produced by.** `Figure_1/cd_recombination_breakpoints/scripts/01_annotate_recombination_dates.py`.

### 1.6 Sequence / structural feature profiles (Figures 2, S3)

* Per-position means are computed over the per-sequence RNA folding predictions:
`Mean_GC`, `Mean_Energy` (ΔG, LinearFold) and the unpaired fraction
`Dot_Proportion`; the averaged paired nucleotide proportion is
`Proppaired = 1 − Dot_Proportion`. Nucleotide diversity π is averaged in the
same way.
* **Windowing.** All features are averaged on a common sliding window of
**100 nt with a 10-nt step** (2,970 windows), so every feature and the
recombination intensity are aligned on one grid.
* **Region classification.** Windows are labelled **Cold** (bottom 10% of
recombination intensity), **Warm** (middle) and **Hot** (top 10%).
* **Statistics.** Two-sided **Mann-Whitney** tests compare Cold with Hot for each
feature (Figure 2b); **Spearman rank** correlations (global and stratified by
region type) relate each feature to the recombination intensity (Figure S3),
with Benjamini-Hochberg correction in both cases.
* **Produced by.** `Figure_2/scripts/01`–`03`, `Figure_S3/scripts/01`.

### 1.7 Spatial-hotspot permutation test (Figure S2h–j)

* For each genomic region, a one-sided Monte-Carlo permutation test
(`perm::permTS`, `exact.mc`, **N = 5,000** resamples, seed 42) asks whether the
mean recombination intensity inside the region exceeds the background formed by
all positions outside it; p-values are adjusted with Benjamini-Hochberg.
* **Produced by.** `Figure_S2/permutation_hotspots/scripts/01_permutation_test.R`.

### 1.8 Phylogenetic diversification (Figure S1)

* The Nextstrain "ncov open / global / all-time" v2 tree JSON (subsampled global
build, 2,938 tips) is traversed recursively and rendered as a rectangular
cladogram; tips are coloured by Nextstrain `clade_membership` and the
recombinant / variant-of-interest tip set is extracted from the matching
metadata.
* **Produced by.** `Figure_S1/scripts/plot_variant_phylogeny.py`.

### 1.9 Structural mapping (Figure 1f–g)

* The per-residue recombination intensity is written into the **B-factor column**
of the PDB 6VXX spike trimer, so the intensity is rendered as a colour ramp
(grey → red) in PyMOL/ChimeraX; the hotspot windows (residues 444–455 and
681–702) and the RBD (438–508) are zoomed for the panel insets.
* **Produced by.** `Figure_1/fg_structural_mapping/map_intensity_to_spike.py`.

\---

## 2\. Data generated in this study

Every table a script writes goes to `data/derived/`
(intermediates), never into `data/raw/` or `data/reference/`.


**Inputs distributed with the code** (`data/raw/`, `data/reference/`) — these are
the tables the figures are drawn from:

|Dataset|Content|Produced by|
|-|-|-|
|`data/reference/recomb_intensity_per_base.csv`|uniformed recombination intensity at 1-nt resolution across the genome, positions 0–29,902|`config.R::uniformed_intensity()` / originally `02_breakpoint_intervals.R` (see §1.3)|
|`data/raw/descendants_node_dates.tsv`|node → earliest descendant collection date, with descendant counts and date quality|`Figure_1/cd_recombination_breakpoints/scripts/01_annotate_recombination_dates.py`|
|`data/raw/hcov-19_proportion.csv`|variant proportion per collection date and location (recombinant groups)|`Figure_1/a_variant_proportion/scripts/02_aggregate_recombvar_proportion.R`|
|`data/raw/hcov-19_VOIproportion.csv`|variant proportion per collection date and location (all variants of interest)|`Figure_1/a_variant_proportion/scripts/03_aggregate_allvar_proportion.R`|
|`data/derived/SARS2_Results.xlsx`|windowed feature matrix: position, GC%, ΔG, Proppaired, π, recombination intensity, sample count, TRS (2,970 windows); drives Figure 2 and Figure S3|`Figure_2/scripts/01`–`02` (+ upstream folding)|

**Intermediates** (`data/derived/`, rebuilt on demand):

|Dataset|Content|Produced by|
|-|-|-|
|`recombination_annotated_with_dates.tsv`|every recombination row + earliest descendant collection date|`01_annotate_recombination_dates.py`|
|`recombination_events_summary.tsv`|one row per recombination node: earliest date, descendant counts, date quality, best `recomb_parsimony`|`01_annotate_recombination_dates.py`|
|`SARS2_summary_results.csv`|per-position Mean\_GC, Mean\_Energy, Dot\_Proportion, Sample\_Count|`Figure_2/scripts/01_feature_summary.R`|
|`SARS2_Feature_Windows.csv`|GC%, ΔG, Proppaired, π and uniformed recombination intensity on the shared 100-nt / 10-nt grid|`Figure_2/scripts/02_feature_windows.R`|

## 3\. Third-party data used (not generated here)

|Data|Source / accession|
|-|-|
|UCSC SARS-CoV-2 public variant set, `public-2025-12-31`|UCSC (used upstream of the recombination calling)|
|Nextstrain "ncov open / global / all-time" tree|Nextstrain|
|PDB 6VXX spike ectodomain trimer|[https://doi.org/10.2210/pdb6VXX/pdb](https://doi.org/10.2210/pdb6VXX/pdb)|
|Cryo-EM maps and models: EMD-69781 / PDB 24RM, EMD-69764 / PDB 24QL, EMD-69824 / PDB 24UF, EMD-69779 / PDB 24RG, EMD-69823 / PDB 24UE, EMD-69778 / PDB 24RE, EMD-69822 / PDB 24UD, EMD-69765 / PDB 24QM|EMDB / PDB (deposited in this study)|
|Pango lineage designations and lineage-proposal issue pages|https://github.com/sars-cov-2-variants/lineage-proposals|



