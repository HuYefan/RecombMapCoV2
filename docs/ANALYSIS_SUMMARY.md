# Analysis summary — algorithms and data generated

This document summarises **what the code computes** (the algorithms) and **what
data it produces** (the generated datasets), figure by figure. It is written to
be usable directly when filling in a Data Availability statement.

Repository: **https://github.com/HuYefan/RecombMapCoV2** (code `RecombMapCoV2`)

\---

## 1\. Algorithms

### 1.1 Variant prevalence (Figure 1a)

* **Lineage grouping.** Each Pango lineage is reduced to its prefix (text before
the first `.`), `(consensus call)` suffixes are stripped, and the prefix is
mapped to a recombinant group through the `Sublineage` field of
`RecombVar\\\_ofInterest.xlsx` / `Variants\\\_of\\\_Interest.xlsx` (a value such as
`XDV|NB` expands to several prefixes). Unmatched lineages become `Others`.
* **Proportions.** Counts are aggregated per collection date and per location
(country), and a global `World` aggregate is appended
(`Proportion = count / Total`).
* **Smoothing.** A **7-day centred rolling sum** is applied to both numerator and
denominator before the ratio is taken
(`Smoothed\\\_Proportion = Smoothed\\\_count / Smoothed\\\_Total`), which avoids gaps
and unstable ratios on low-count days.
* **Produced by.** `Figure\\\_1/a\\\_variant\\\_proportion/scripts/01`–`05`.

### 1.2 Recombination-pathway delineation (Figure 1b)

* For each recombinant lineage the lineage-defining sequences and the two
parental references are MAFFT-aligned (`--auto`) and drawn in snipit
**recombination mode** (`--recombi-mode --recombi-references`), which assigns
every alignment column to the parental lineage it most resembles and infers the
switch points.
* **Produced by.** `Figure\\\_1/b\\\_recombination\\\_pathways/recombination\\\_pathways.py`.

### 1.3 Uniformed recombination intensity (Figures 1c–1e, 2a, S2)

* **Event filter.** Only recombination nodes with `recomb\\\_parsimony >= 7` are
used (1,273 nodes in RIPPLES output; 55,498 breakpoint intervals).
* **Interval weighting.** A recombination event localises each breakpoint to an
interval between two recombination-informative SNPs, so the event is spread
**uniformly** over that interval: an interval of length *L* contributes weight
**1/L** to each of its *L* bases. The per-base **uniformed recombination
intensity** is the sum over all intervals covering a base.
* **Implementation.** A difference array with one cumulative sum, with the
weights aggregated per index (`rowsum`) — an O(n) equivalent of the per-base
loop. In R a naive `delta\\\[idx] <- delta\\\[idx] + weight` is not an accumulating
update when the subscript repeats, so the weights must be summed per index
first. The track the submitted figures were drawn from ships as
`data/reference/recomb\\\_intensity\\\_per\\\_base.csv` and is what every panel reads
through `config.R::load\\\_recomb\\\_intensity()`; the re-derivation is kept in
`config.R::uniformed\\\_intensity()` and covered by
`tests/test\\\_uniformed\\\_intensity.R`.
* **Peak regions.** nt ≈ 284–302 and the RBM window nt 22,895–22,927 (S:445–455).
* **Amino-acid resolution.** The 1-nt intensity restricted to the spike CDS
(21,563–25,384) is summed over each codon, giving a per-residue intensity for
S1–S1274.
* **Threshold.** mean + 3 SD of the genome-wide per-base intensity (dashed line
in the figures).
* **Peak regions (both tracks).** nt ≈ 284–302 and the RBM window
nt 22,895–22,927 (S:445–455).
* **Regression test.** `tests/test\\\_uniformed\\\_intensity.R` (8 checks, including a
comparison against the straightforward per-base reference loop on the real
event table).
* **Produced by.** `Figure\\\_1/cd\\\_recombination\\\_breakpoints/scripts/02\\\_breakpoint\\\_intervals.R`,
`Figure\\\_1/e\\\_recombination\\\_intensity/scripts/01\\\_plot\\\_spike\\\_intensity.R`,
`Figure\\\_S2/uniformed\\\_intensity/scripts/01\\\_uniformed\\\_intensity.R`.

### 1.4 Per-base midpoint histogram (Figure 1d)

* Because an event can only be localised to an interval, the **midpoint** of each
breakpoint interval is used as the point estimate of the breakpoint position;
the histogram uses a **bin width of 1 nt**.
* **Produced by.** `Figure\\\_1/cd\\\_recombination\\\_breakpoints/scripts/02\\\_breakpoint\\\_intervals.R`.

### 1.5 Event dating (Figure S2, supporting panels)

* Each recombination node is dated by the **earliest collection date among its
descendant sequences** (an upper bound on the true recombination time). Only
complete `YYYY-MM-DD` dates are used for the primary estimate; nodes without
one fall back to partial dates (`YYYY-MM-01`, `YYYY-01-01`) and are flagged
`partial\\\_only`.
* **Produced by.** `Figure\\\_1/cd\\\_recombination\\\_breakpoints/scripts/01\\\_annotate\\\_recombination\\\_dates.py`.

### 1.6 Sequence / structural feature profiles (Figures 2, S3)

* Per-position means are computed over the per-sequence RNA folding predictions:
`Mean\\\_GC`, `Mean\\\_Energy` (ΔG, LinearFold) and the unpaired fraction
`Dot\\\_Proportion`; the averaged paired nucleotide proportion is
`Proppaired = 1 − Dot\\\_Proportion`. Nucleotide diversity π is averaged in the
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
* **Produced by.** `Figure\\\_2/scripts/01`–`03`, `Figure\\\_S3/scripts/01`.

### 1.7 Spatial-hotspot permutation test (Figure S2h–j)

* For each genomic region, a one-sided Monte-Carlo permutation test
(`perm::permTS`, `exact.mc`, **N = 5,000** resamples, seed 42) asks whether the
mean recombination intensity inside the region exceeds the background formed by
all positions outside it; p-values are adjusted with Benjamini-Hochberg.
* **Produced by.** `Figure\\\_S2/permutation\\\_hotspots/scripts/01\\\_permutation\\\_test.R`.

### 1.8 Phylogenetic diversification (Figure S1)

* The Nextstrain "ncov open / global / all-time" v2 tree JSON (subsampled global
build, 2,938 tips) is traversed recursively and rendered as a rectangular
cladogram; tips are coloured by Nextstrain `clade\\\_membership` and the
recombinant / variant-of-interest tip set is extracted from the matching
metadata.
* **Produced by.** `Figure\\\_S1/scripts/plot\\\_variant\\\_phylogeny.py`.

### 1.9 Structural mapping (Figure 1f–g)

* The per-residue recombination intensity is written into the **B-factor column**
of the PDB 6VXX spike trimer, so the intensity is rendered as a colour ramp
(grey → red) in PyMOL/ChimeraX; the hotspot windows (residues 444–455 and
681–702) and the RBD (438–508) are zoomed for the panel insets.
* **Produced by.** `Figure\\\_1/fg\\\_structural\\\_mapping/map\\\_intensity\\\_to\\\_spike.py`.

\---

## 2\. Data generated in this study

Every table a script writes goes to `figures/` (figure outputs) or `data/derived/`
(intermediates), never into `data/raw/` or `data/reference/`.

**Tables written next to the figures** (`figures/main/`, `figures/supplementary/`):

|Generated dataset|Content|Rows|Produced by|
|-|-|-|-|
|`figures/main/Figure1e\\\_spike\\\_intensity\\\_aa.csv`|uniformed recombination intensity per spike residue (S1–S1274)|1,274|`Figure\\\_1/e\\\_recombination\\\_intensity/scripts/01\\\_plot\\\_spike\\\_intensity.R`|
|`figures/main/Figure2b\\\_cold\\\_hot\\\_statistics.csv`|Mann-Whitney results (W, p, BH-adjusted p) per feature, Cold vs. Hot|5|`Figure\\\_2/scripts/03\\\_plot\\\_genome\\\_features.R`|
|`figures/supplementary/FigureS3\\\_pairwise\\\_correlations.csv`|Spearman ρ, p and BH-adjusted p for every feature × intensity, globally and per region type|20|`Figure\\\_S3/scripts/01\\\_vcm\\\_correlation\\\_matrix.R`|
|`figures/supplementary/FigureS2h\\\_permutation\\\_results.tsv`|permutation-test Z-scores, raw and BH-adjusted p per genomic region|37|`Figure\\\_S2/permutation\\\_hotspots/scripts/01\\\_permutation\\\_test.R`|
|`figures/supplementary/FigureS2j\\\_uniformed\\\_intensity\\\_spike\\\_aa.csv`|per-residue uniformed intensity (Figure S2j panel values)|1,274|`Figure\\\_S2/uniformed\\\_intensity/scripts/01\\\_uniformed\\\_intensity.R`|
|`figures/supplementary/FigureS1b\\\_variant\\\_tip\\\_counts.csv`|number of Nextstrain tips per variant of interest|8|`Figure\\\_S1/scripts/plot\\\_variant\\\_phylogeny.py`|
|`figures/main/Figure1fg\\\_pymol.pml`, `Figure1f\\\_6vxx\\\_intensity\\\_bfactor.pdb`|PyMOL script and PDB 6VXX with the per-residue intensity in the B-factor column|—|`Figure\\\_1/fg\\\_structural\\\_mapping/map\\\_intensity\\\_to\\\_spike.py`|

**Inputs distributed with the code** (`data/raw/`, `data/reference/`) — these are
the tables the figures are drawn from:

|Dataset|Content|Produced by|
|-|-|-|
|`data/reference/recomb\\\_intensity\\\_per\\\_base.csv`|uniformed recombination intensity at 1-nt resolution across the genome, positions 0–29,902|`config.R::uniformed\\\_intensity()` / originally `02\\\_breakpoint\\\_intervals.R` (see §1.3)|
|`data/raw/descendants\\\_node\\\_dates.tsv`|node → earliest descendant collection date, with descendant counts and date quality|`Figure\\\_1/cd\\\_recombination\\\_breakpoints/scripts/01\\\_annotate\\\_recombination\\\_dates.py`|
|`data/raw/hcov-19\\\_proportion.csv`|variant proportion per collection date and location (recombinant groups)|`Figure\\\_1/a\\\_variant\\\_proportion/scripts/02\\\_aggregate\\\_recombvar\\\_proportion.R`|
|`data/raw/hcov-19\\\_VOIproportion.csv`|variant proportion per collection date and location (all variants of interest)|`Figure\\\_1/a\\\_variant\\\_proportion/scripts/03\\\_aggregate\\\_allvar\\\_proportion.R`|
|`data/derived/SARS2\\\_Results.xlsx`|windowed feature matrix: position, GC%, ΔG, Proppaired, π, recombination intensity, sample count, TRS (2,970 windows); drives Figure 2 and Figure S3|`Figure\\\_2/scripts/01`–`02` (+ upstream folding)|

**Intermediates** (`data/derived/`, rebuilt on demand):

|Dataset|Content|Produced by|
|-|-|-|
|`recombination\\\_annotated\\\_with\\\_dates.tsv`|every recombination row + earliest descendant collection date|`01\\\_annotate\\\_recombination\\\_dates.py`|
|`recombination\\\_events\\\_summary.tsv`|one row per recombination node: earliest date, descendant counts, date quality, best `recomb\\\_parsimony`|`01\\\_annotate\\\_recombination\\\_dates.py`|
|`SARS2\\\_summary\\\_results.csv`|per-position Mean\_GC, Mean\_Energy, Dot\_Proportion, Sample\_Count|`Figure\\\_2/scripts/01\\\_feature\\\_summary.R`|
|`SARS2\\\_Feature\\\_Windows.csv`|GC%, ΔG, Proppaired, π and uniformed recombination intensity on the shared 100-nt / 10-nt grid|`Figure\\\_2/scripts/02\\\_feature\\\_windows.R`|

## 3\. Third-party data used (not generated here)

|Data|Source / accession|
|-|-|
|UCSC SARS-CoV-2 public variant set, `public-2025-12-31`|UCSC (used upstream of the recombination calling)|
|Nextstrain "ncov open / global / all-time" tree|Nextstrain|
|PDB 6VXX spike ectodomain trimer|[https://doi.org/10.2210/pdb6VXX/pdb](https://doi.org/10.2210/pdb6VXX/pdb)|
|Cryo-EM maps and models: EMD-69781 / PDB 24RM, EMD-69764 / PDB 24QL, EMD-69824 / PDB 24UF, EMD-69779 / PDB 24RG, EMD-69823 / PDB 24UE, EMD-69778 / PDB 24RE, EMD-69822 / PDB 24UD, EMD-69765 / PDB 24QM|EMDB / PDB (deposited in this study)|
|Pango lineage designations and lineage-proposal issue pages|https://github.com/sars-cov-2-variants/lineage-proposals|



