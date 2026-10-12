# Figure 1 — global prevalence and breakpoint analysis of SARS-CoV-2 recombinants

|Panel|Content|Folder|
|-|-|-|
|a|Area charts of the relative proportion of dominant recombinant variants (XBB, XEC, XFG, XDV) against non-recombinant variants, globally and for Asia / Europe / North America|[`a\\\\\\\_variant\\\\\\\_proportion`](a_variant_proportion)|
|b|Recombination pathways and the spike mutations acquired by XBB, XDV, XEC and XFG from their parental lineages|[`b\\\\\\\_recombination\\\\\\\_pathways`](b_recombination_pathways)|
|c|Genome-wide mapping of the 2,991 recombination breakpoint intervals|[`cd\\\\\\\_recombination\\\\\\\_breakpoints`](cd_recombination_breakpoints)|
|d|Midpoint histogram of the breakpoint intervals (bin = 1 nt)|[`cd\\\\\\\_recombination\\\\\\\_breakpoints`](cd_recombination_breakpoints)|
|e|High-resolution uniformed recombination intensity|[`e\\\\\\\_recombination\\\\\\\_intensity`](e_recombination_intensity)|
|f|Uniformed recombination intensity mapped onto the WT spike structure|[`fg\\\\\\\_structural\\\\\\\_mapping`](fg_structural_mapping)|
|g|Zoomed views of the hotspot regions in the spike protein|[`fg\\\\\\\_structural\\\\\\\_mapping`](fg_structural_mapping)|

## Run order

```bash
# a) variant proportion, World / Asia / Europe / North America
Rscript Figure\\\\\\\_1/a\\\\\\\_variant\\\\\\\_proportion/scripts/01\\\\\\\_prepare\\\\\\\_metadata.R            # needs the GISAID dump
Rscript Figure\\\\\\\_1/a\\\\\\\_variant\\\\\\\_proportion/scripts/02\\\\\\\_aggregate\\\\\\\_recombvar\\\\\\\_proportion.R
Rscript Figure\\\\\\\_1/a\\\\\\\_variant\\\\\\\_proportion/scripts/03\\\\\\\_aggregate\\\\\\\_allvar\\\\\\\_proportion.R
Rscript Figure\\\\\\\_1/a\\\\\\\_variant\\\\\\\_proportion/scripts/04\\\\\\\_plot\\\\\\\_recombvar\\\\\\\_proportion.R   # -> Figure 1a
Rscript Figure\\\\\\\_1/a\\\\\\\_variant\\\\\\\_proportion/scripts/05\\\\\\\_plot\\\\\\\_allvar\\\\\\\_proportion.R

# c-e) breakpoints and uniformed intensity
python Figure\\\\\\\_1/cd\\\\\\\_recombination\\\\\\\_breakpoints/scripts/01\\\\\\\_annotate\\\\\\\_recombination\\\\\\\_dates.py
Rscript Figure\\\\\\\_1/cd\\\\\\\_recombination\\\\\\\_breakpoints/scripts/02\\\\\\\_breakpoint\\\\\\\_intervals.R  # -> Figure 1c, 1d, 1e
Rscript Figure\\\\\\\_1/e\\\\\\\_recombination\\\\\\\_intensity/scripts/01\\\\\\\_plot\\\\\\\_spike\\\\\\\_intensity.R     # -> Figure 1e (per AA)

# f-g) structural mapping
python Figure\\\\\\\_1/fg\\\\\\\_structural\\\\\\\_mapping/map\\\\\\\_intensity\\\\\\\_to\\\\\\\_spike.py --pymol

# b) recombination pathway diagrams (needs mafft + snipit)
python Figure\\\\\\\_1/b\\\\\\\_recombination\\\\\\\_pathways/recombination\\\\\\\_pathways.py --list
python Figure\\\\\\\_1/b\\\\\\\_recombination\\\\\\\_pathways/recombination\\\\\\\_pathways.py --fasta-dir fasta
```

Steps 01–03 of panel a only need to be re-run when the underlying GISAID dump
changes; the aggregated tables are distributed and steps 04/05 run directly.

## Inputs

|File|Distributed|Description|
|-|-|-|
|`../data/raw/hcov-19\\\\\\\_proportion.csv`|yes (1.8 MB)|variant proportions per date / region (recombinant panel)|
|`../data/raw/hcov-19\\\\\\\_VOIproportion.csv`|yes (3.0 MB)|variant proportions per date / region (all variants of interest)|
|`../data/raw/recombination.tsv`|yes (9.0 MB)|RIPPLES recombination events, 1,273 nodes|
|`../data/raw/descendants\\\\\\\_node\\\\\\\_dates.tsv`|yes (47 KB)|earliest descendant collection date per node|
|`../data/reference/recomb\\\\\\\_intensity\\\\\\\_per\\\\\\\_base.csv`|yes (541 KB)|per-base uniformed recombination intensity (panels c–e)|
|`../data/reference/RecombVar\\\\\\\_ofInterest.xlsx`|yes|lineage → sublineage mapping (recombinant groups)|
|`../data/reference/Variants\\\\\\\_of\\\\\\\_Interest.xlsx`|yes|lineage → sublineage mapping (all variants)|
|`../data/reference/SARS-CoV-2\\\\\\\_ORF.xlsx`, `SARS-CoV-2\\\\\\\_Spike.xlsx`|yes|ORF / spike sub-domain coordinates|
|`../data/reference/Recombinant\\\\\\\_Variants.xlsx`|yes|designated recombinants and their breakpoints|
|`../data/reference/6vxx.pdb`|yes (4.7 MB)|spike ectodomain trimer (PDB 6VXX)|



