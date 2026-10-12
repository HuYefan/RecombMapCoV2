#!/usr/bin/env python3
# =============================================================================
#  map_intensity_to_spike.py
#  Figure 1f-g : structural mapping of the uniformed recombination intensity on
#  the SARS-CoV-2 spike trimer, plus zoomed views of the hotspot regions.
# =============================================================================
#
#  What it does
#  ------------
#    1. reads figures/main/Figure1e_spike_intensity_aa.csv (per-residue uniformed
#       recombination intensity, produced by Figure_1/e_.../01_plot_spike_intensity.R)
#    2. writes a per-residue B-factor column into a copy of the spike PDB, so the
#       intensity can be shown as a colour ramp in any molecular viewer
#    3. optionally drives PyMOL to render the overview (panel f) and the zoomed
#       hotspot views (panel g)
#
#  Input
#  -----
#    data/reference/6vxx.pdb      spike ectodomain trimer (closed, one-RBD-up not
#                                 modelled here - the manuscript renders the WT
#                                 spike from 6VXX)
#    figures/main/Figure1e_spike_intensity_aa.csv
#
#  Output
#  ------
#    figures/main/Figure1f_6vxx_intensity_bfactor.pdb
#    figures/main/Figure1f_spike_intensity.png
#    figures/main/Figure1g_spike_hotspots.png      (requires PyMOL)
#
#  Usage
#  -----
#    python Figure_1/fg_structural_mapping/map_intensity_to_spike.py
#    python Figure_1/fg_structural_mapping/map_intensity_to_spike.py --pymol
# =============================================================================

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config import DATA, DATA_REF, DERIVED, FIG_MAIN  # noqa: E402

# SARS-CoV-2 spike residues that carry a recombination breakpoint hotspot in the
# manuscript (RBM 445-455 and the S1/S2 region); used for the zoomed panels.
HOTSPOT_RESIDUES = [
    (444, 455, "RBM hotspot (S:445-455)"),
    (681, 702, "S1/S2 cleavage site"),
]
COLOUR_LOW = "grey80"
COLOUR_HIGH = "red"


def read_intensity(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            values[int(float(row["aa_pos"]))] = float(row["intensity"])
    return values


def write_bfactors(pdb_in: Path, pdb_out: Path, intensity: dict[int, float]) -> None:
    """Copy pdb_in to pdb_out, replacing the B-factor of every spike residue
    with its uniformed recombination intensity (0 where no data)."""
    scale = max(intensity.values()) if intensity else 1.0
    n_written = 0
    with open(pdb_in, encoding="utf-8") as fi, open(pdb_out, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    fo.write(line)
                    continue
                val = intensity.get(resseq, 0.0)
                # B-factor occupies columns 61-66 (0-based 60:66)
                line = f"{line[:60]}{val:6.2f}{line[66:]}"
                n_written += 1
            fo.write(line)
    print(f"wrote {pdb_out}  ({n_written} atoms; scale 1.0 = {scale:.3f})")


PYMOL_HEADER = """
# ---- Figure 1f: overview of the spike trimer --------------------------------
load {pdb}, spike
hide everything
bg_color white
set ray_shadows, 0
set antialias, 2

# per-residue intensity -> spectrum, grey where there is no data
spectrum b, {low} {high}, spike, minimum=0, maximum={vmax}
cartoon automatic, spike
show cartoon, spike

# RBD highlighted for orientation
select rbd, spike and resi 438-508
color marine, rbd

orient spike
set ray_opaque_background, 0
png {out_f}, width=1800, dpi=300, ray=1

# ---- Figure 1g: zoomed hotspot views ---------------------------------------
"""

PYMOL_HOTSPOT = """
# {label}
select hot, spike and resi {a}-{b}
zoom hot, buffer=3
png {out_g_prefix}_{a}_{b}.png, width=1200, dpi=300, ray=1
"""


def render_with_pymol(pdb: Path, intensity: dict[int, float], hotspots: list[tuple[int, int, str]]) -> None:
    exe = shutil.which("pymol")
    if exe is None:
        print("[skip] PyMOL not found on PATH; the B-factor PDB can be rendered "
              "manually (colour by B-factor, spectrum grey->red).")
        return

    out_f = FIG_MAIN / "Figure1f_spike_intensity.png"
    out_g_prefix = FIG_MAIN / "Figure1g_spike_hotspot"

    script = PYMOL_HEADER.format(
        pdb=str(pdb).replace("\\", "/"),
        low=COLOUR_LOW,
        high=COLOUR_HIGH,
        vmax=max(intensity.values()) if intensity else 1.0,
        out_f=str(out_f).replace("\\", "/"),
    )
    for a, b, label in hotspots:
        script += PYMOL_HOTSPOT.format(
            label=label, a=a, b=b,
            out_g_prefix=str(out_g_prefix).replace("\\", "/"),
        )

    script_file = FIG_MAIN / "Figure1fg_pymol.pml"
    script_file.write_text(script, encoding="utf-8")
    print(f"wrote {script_file}")
    subprocess.run([exe, "-cq", str(script_file)], check=True)
    print(f"rendered -> {out_f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pymol", action="store_true",
                    help="also render the panels with PyMOL (must be on PATH)")
    ap.add_argument("--pdb", default=str(DATA_REF / "6vxx.pdb"),
                    help="input spike PDB (default: data/reference/6vxx.pdb)")
    args = ap.parse_args()

    pdb_in = Path(args.pdb)
    if not pdb_in.exists():
        sys.exit(f"missing {pdb_in}\nSee data/README.md for the 6VXX structure.")

    intensity_file = FIG_MAIN / "Figure1e_spike_intensity_aa.csv"
    if not intensity_file.exists():
        sys.exit(f"missing {intensity_file}\nRun Figure_1/e_recombination_intensity/"
                 f"scripts/01_plot_spike_intensity.R first.")

    intensity = read_intensity(intensity_file)
    print(f"intensity values for {len(intensity)} spike residues; "
          f"max = {max(intensity.values()):.3f}")

    pdb_out = FIG_MAIN / "Figure1f_6vxx_intensity_bfactor.pdb"
    write_bfactors(pdb_in, pdb_out, intensity)

    if args.pymol:
        render_with_pymol(pdb_out, intensity, HOTSPOT_RESIDUES)
    else:
        print("\nNext: open figures/main/Figure1f_6vxx_intensity_bfactor.pdb and colour by "
              "B-factor (grey -> red) for Figure 1f, then zoom the residues "
              + ", ".join(f"{a}-{b}" for a, b, _ in HOTSPOT_RESIDUES) + " for Figure 1g.")


if __name__ == "__main__":
    main()
