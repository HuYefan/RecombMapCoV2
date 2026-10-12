# Figure 1a — variant proportion over time

Area charts showing the relative proportion of the dominant recombinant
variants (XBB, XEC, XFG, XDV) compared with non-recombinant variants (light
grey), globally and for Asia, Europe and North America. A companion panel shows
the full variant-of-interest set (Alpha → XFG).

## Scripts (run in order, from the repository root)

|Step|Script|Input|Output|
|-|-|-|-|
|1|`scripts/01\_prepare\_metadata.R`|raw metadata dump (not distributed)|`../data/raw/hcov-19\_metadata\_filtered.csv`, `../data/raw/hcov-19\_AA\_Substitutions.csv`|
|2|`scripts/02\_aggregate\_recombvar\_proportion.R`|step 1 + `RecombVar\_ofInterest.xlsx`|`../data/raw/hcov-19\_proportion.csv`|
|3|`scripts/03\_aggregate\_allvar\_proportion.R`|step 1 + `Variants\_of\_Interest.xlsx`|`../data/raw/hcov-19\_VOIproportion.csv`|
|4|`scripts/04\_plot\_recombvar\_proportion.R`|step 2|**Figure 1a** (area + line)|
|5|`scripts/05\_plot\_allvar\_proportion.R`|step 3|all-variant area chart|

Steps 2 and 3 outputs are distributed with the repository, so the plotting
scripts (4 and 5) run without the multi-gigabyte GISAID dump.

## Method notes

* A lineage is assigned to a recombinant group through the `Sublineage` column of
`RecombVar\_ofInterest.xlsx` / `Variants\_of\_Interest.xlsx`: a value such as
`"XDV|NB"` expands into the prefixes `XDV` and `NB`, matched against the part of
the Pango lineage before the first dot. Anything unmatched is `Others`.
* Proportions are computed per collection date **and** per location, and a global
`World` row set is appended.
* The plotted curves are a **7-day centred rolling mean** of the daily counts, so
the numerator and denominator are smoothed before the ratio is taken
(`Smoothed\_count / Smoothed\_Total`).
* *(Consensus call)* suffixes are stripped from the Pango lineage before matching.

## Outputs

* `../figures/main/Figure1a\_recombvar\_proportion\_area.pdf`
* `../figures/main/Figure1a\_recombvar\_proportion\_line.pdf`
* `../figures/main/Figure1a\_allvar\_proportion\_area.pdf`

