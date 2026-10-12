# Figure 1f–g — structural mapping of recombination intensity on the spike

**Panel f.** Uniformed recombination intensity mapped onto the SARS-CoV-2
wild-type spike structure.

**Panel g.** Zoomed views of the highlighted hotspot regions within the spike
protein (red and blue denote the uniformed recombination intensity of spike
regions).

## Script

| Script | Purpose |
|---|---|
| `map_intensity_to_spike.py` | write the per-residue intensity into the B-factor column of the spike PDB and (optionally) render panel f and the zoomed panel g with PyMOL |

```bash
# 1. produce the per-residue intensity (once)
Rscript Figure_1/e_recombination_intensity/scripts/01_plot_spike_intensity.R

# 2. write the B-factor PDB
python Figure_1/fg_structural_mapping/map_intensity_to_spike.py

# 3. render with PyMOL (needs `pymol` on PATH)
python Figure_1/fg_structural_mapping/map_intensity_to_spike.py --pymol
```

## Inputs

| File | Description |
|---|---|
| `../../data/reference/6vxx.pdb` | **distributed** — spike ectodomain trimer, PDB [6VXX](https://www.rcsb.org/structure/6VXX) |
| `../../figures/main/Figure1e_spike_intensity_aa.csv` | per-residue uniformed intensity (written by the Figure 1e step) |

## Outputs

* `../../figures/main/Figure1f_6vxx_intensity_bfactor.pdb` — the structure with the
  intensity in the B-factor column (colour by B-factor, grey → red)
* `../../figures/main/Figure1fg_pymol.pml` — the generated PyMOL script
* `../../figures/main/Figure1f_spike_intensity.png`
* `../../figures/main/Figure1g_spike_hotspot_<start>_<end>.png`

## Hotspot residues

Defined by `HOTSPOT_RESIDUES` in the script:

| Range | Region |
|---|---|
| 444–455 | RBM hotspot (S:445–455), the hottest recombination breakpoint window |
| 681–702 | S1/S2 cleavage site |

If PyMOL is not installed, open
`../../figures/main/Figure1f_6vxx_intensity_bfactor.pdb` in ChimeraX/PyMOL, colour by
B-factor with a grey→red ramp, and orient on residues 438–508 (RBD) for
panel f.
