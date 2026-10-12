#!/usr/bin/env python3
# =============================================================================
#  recombination_pathways.py
#  Figure 1b - recombination-pathway diagrams (snipit plots) for XBB, XDE, XDV,
#  XEC and XFG, plus the corresponding reference ("--Ref") panels.
# =============================================================================
#
#  Method
#  ------
#  For each recombinant lineage, the designated lineage-defining sequences were
#  downloaded from the Pango lineage-proposal pages (the exact accessions and
#  the source issue for every lineage are listed in LINEAGES below), aligned
#  with MAFFT and drawn in "recombination mode" with snipit, using the two
#  parental lineages as the recombination references.  The resulting PDFs are
#  the panels assembled as Figure 1b.
#
#  Requirements
#  ------------
#      conda install -c bioconda mafft snipit
#      pip install biopython
#
#  Usage
#  -----
#      # 1. put the FASTA (recombinant + both parents + Ref) in a folder
#      # 2. run
#      python Figure_1/b_recombination_pathways/recombination_pathways.py \
#          --fasta-dir fasta --outdir figures/main/Figure1b
#
#  The original notebook ran on Google Colab with interactive file uploads; this
#  script is the same analysis with the file handling made non-interactive.
# =============================================================================

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Provenance for every lineage drawn in Figure 1b
# ---------------------------------------------------------------------------
LINEAGES = {
    "XBB": {
        "reference_page": "https://github.com/cov-lineages/pango-designation/issues/1058",
        "parents": ["BM.1.1.1", "BJ.1"],
        "sequences": {
            "BJ.1": "hCoV-19/India/MH-GMCA-ICMR-INSACOG-C635137/2022|EPI_ISL_14923255|2022-08-29",
            "BM.1.1.1": "hCoV-19/India/MH-GMCA-ICMR-INSACOG-C635652/2022|EPI_ISL_14923234|2022-09-02",
            "XBB": "hCoV-19/Singapore/11487/2022|EPI_ISL_14917761|2022-09-06",
        },
    },
    "XDE": {
        "reference_page": "https://github.com/sars-cov-2-variants/lineage-proposals/issues/855",
        "parents": ["FL.13.4", "GW.5.1"],
        "sequences": {
            "FL.13.4": "EPI_ISL_17690187",
            "GW.5.1": "EPI_ISL_18112250",
            "XDE": "EPI_ISL_18294822",
        },
    },
    "XDV": {
        "reference_page": "https://github.com/cov-lineages/pango-designation/issues/2402",
        "parents": ["XDE", "JN.1"],
        "sequences": {
            "JN.1": "EPI_ISL_18550520",
            "XDE": "EPI_ISL_18389828",
            "XDV": "EPI_ISL_18576372",
        },
    },
    "XEC": {
        "reference_page": ("https://github.com/sars-cov-2-variants/lineage-proposals/issues/1827 / "
                           "https://github.com/cov-lineages/pango-designation/issues/2717"),
        "parents": ["KS.1.1", "KP.3.3"],
        "sequences": {
            "KS.1.1": "EPI_ISL_19265606",
            "KP.3.3": "EPI_ISL_19280812",
            "XEC": "EPI_ISL_19283881",
        },
    },
    "XFG": {
        "reference_page": "https://github.com/sars-cov-2-variants/lineage-proposals/issues/2449",
        "parents": ["LF.7", "NW.1.1.1"],
        "sequences": {
            "LF.7": "EPI_ISL_19698555",
            "NW.1.1.1": "EPI_ISL_19810352",
            "XFG": "EPI_ISL_19724626",
        },
    },
}

COLOUR_PALETTE = "classic_extended"


def which_or_die(program: str) -> str:
    path = shutil.which(program)
    if path is None:
        sys.exit(f"'{program}' was not found on PATH. Install with: "
                 f"conda install -c bioconda mafft snipit")
    return path


def remove_sequence_by_id(input_file: Path, output_file: Path, seq_id: str) -> bool:
    """Write input_file's records to output_file, dropping record.id == seq_id."""
    kept, removed = [], False
    for record in SeqIO.parse(str(input_file), "fasta"):
        if record.id != seq_id:
            kept.append(record)
        else:
            removed = True
    if removed:
        SeqIO.write(kept, str(output_file), "fasta")
        print(f"  removed '{seq_id}' -> {output_file.name}")
    return removed


def run(cmd: list[str]) -> None:
    print("  $", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_lineage(name: str, spec: dict, fasta_dir: Path, outdir: Path) -> None:
    fasta = fasta_dir / f"{name}.fasta"
    if not fasta.exists():
        print(f"[skip] {name}: {fasta} not found")
        return

    print(f"[{name}] parents = {', '.join(spec['parents'])}")
    aligned = fasta_dir / f"aligned_{fasta.stem}.fasta"
    # MAFFT alignment (written to file, then reused by both snipit calls)
    with open(aligned, "w") as fh:
        subprocess.run(["mafft", "--auto", str(fasta)], check=True,
                       stdout=fh)

    refs = ",".join(spec["parents"])

    # (1) reference panel (keeps the "Ref" sequence)
    run(["snipit", str(aligned), "--reference", "Ref",
         "--recombi-mode", "--recombi-references", refs,
         "--format", "pdf", "--output-file", str(outdir / f"{name}-Ref"),
         "--colour-palette", COLOUR_PALETTE])

    # (2) query panel (Ref removed, reference is the recombinant itself)
    filtered = fasta_dir / f"filtered_{aligned.stem}.fasta"
    if remove_sequence_by_id(aligned, filtered, "Ref") or filtered.exists():
        run(["snipit", str(filtered), "--reference", name,
             "--recombi-mode", "--recombi-references", refs,
             "--format", "pdf", "--output-file", str(outdir / name),
             "--colour-palette", COLOUR_PALETTE])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-dir", default="fasta",
                    help="folder holding <LINEAGE>.fasta (default: fasta)")
    ap.add_argument("--outdir", default="figures/main/Figure1b",
                    help="output folder for the snipit PDFs")
    ap.add_argument("--lineages", nargs="*", default=list(LINEAGES),
                    help="subset of lineages to draw (default: all)")
    ap.add_argument("--list", action="store_true",
                    help="print the expected accessions and exit")
    args = ap.parse_args()

    if args.list:
        for name, spec in LINEAGES.items():
            print(f"\n{name}  ({spec['reference_page']})")
            print(f"  parents: {', '.join(spec['parents'])}")
            for lin, acc in spec["sequences"].items():
                print(f"    {lin:<12} {acc}")
        return

    which_or_die("mafft")
    which_or_die("snipit")

    fasta_dir = Path(args.fasta_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name in args.lineages:
        if name not in LINEAGES:
            sys.exit(f"unknown lineage '{name}'; choose from {list(LINEAGES)}")
        build_lineage(name, LINEAGES[name], fasta_dir, outdir)

    print(f"\ndone - PDFs written to {outdir.resolve()}")


if __name__ == "__main__":
    main()
