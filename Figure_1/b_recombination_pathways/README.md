# Figure 1b — recombination pathways of the major recombinant lineages

Schematic of the recombination pathways and the spike mutations acquired by
XBB, XDV, XEC and XFG from their parental lineages.

## Script

| Script | Purpose |
|---|---|
| `recombination_pathways.py` | MAFFT-align each lineage's FASTA and draw the snipit recombination-mode panels (query + reference) |

```bash
python Figure_1/b_recombination_pathways/recombination_pathways.py --list
python Figure_1/b_recombination_pathways/recombination_pathways.py --fasta-dir fasta --outdir figures/main/Figure1b
```

## Requirements

```bash
conda install -c bioconda mafft snipit
pip install biopython
```

## Lineages and accessions

The script carries the full provenance table (`LINEAGES`) — for every lineage it
records the Pango lineage-proposal issue the sequences came from, the two
parental lineages used as recombination references and the exact GISAID
accessions of the parental and recombinant sequences:

| Lineage | Parents | Source |
|---|---|---|
| XBB | BM.1.1.1, BJ.1 | [pango-designation#1058](https://github.com/cov-lineages/pango-designation/issues/1058) |
| XDE | FL.13.4, GW.5.1 | [lineage-proposals#855](https://github.com/sars-cov-2-variants/lineage-proposals/issues/855) |
| XDV | XDE, JN.1 | [pango-designation#2402](https://github.com/cov-lineages/pango-designation/issues/2402) |
| XEC | KS.1.1, KP.3.3 | [lineage-proposals#1827](https://github.com/sars-cov-2-variants/lineage-proposals/issues/1827), [pango-designation#2717](https://github.com/cov-lineages/pango-designation/issues/2717) |
| XFG | LF.7, NW.1.1.1 | [lineage-proposals#2449](https://github.com/sars-cov-2-variants/lineage-proposals/issues/2449) |

Pre-aligned FASTA files are **not** distributed — they contain GISAID sequence
data and must be re-downloaded through GISAID under its terms of use.

## Input

One FASTA per lineage in `--fasta-dir`, named `<LINEAGE>.fasta`, containing the
recombinant, both parental reference sequences and a sequence named `Ref`.

## Outputs

`figures/main/Figure1b/<LINEAGE>-Ref.pdf` (reference panel) and
`figures/main/Figure1b/<LINEAGE>.pdf` (query panel).

## Note on the original notebook

This script is the non-interactive version of `Results/Recombination_Case/Recombination.ipynb`,
which ran on Google Colab with `files.upload()`.
