# -*- coding: utf-8 -*-
"""
Annotate each recombination event with the time (earliest descendant
collection date) at which the event is estimated to have occurred.

Method (same as the manuscript pipeline):
  For each recombination node, parse the collection date of every
  descendant sequence in descendants.tsv and take the earliest one.
  This is an upper bound on the true recombination time (the event must
  predate the earliest sampled descendant).
  Only complete dates (YYYY-MM-DD) are used for the primary estimate;
  nodes with no complete dates are filled from partial dates
  (YYYY-MM -> YYYY-MM-01, YYYY -> YYYY-01-01).

Inputs:
  data/raw/recombination.tsv   - recombination events (rows = candidate breakpoints)
  data/raw/descendants.tsv     - node_id -> comma-separated strain|accession|date
                             (669 MB; not distributed with the repository)

Outputs:
  data/derived/recombination_annotated_with_dates.tsv - every recombination row + earliest_date
  data/derived/recombination_events_summary.tsv       - one row per unique event (node)
  data/raw/descendants_node_dates.tsv                     - compact node -> date table
                                                        used by the R code

Run from the repository root:
    python Figure_1/cd_recombination_breakpoints/scripts/01_annotate_recombination_dates.py
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

# repository root = three levels up from this file, then source config.py
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from config import DATA, DERIVED, RECOMB_PARSIMONY_MIN  # noqa: E402

RECOMB_FILE = DATA / "recombination.tsv"
DESC_FILE = DATA / "descendants.tsv"
ANNOTATED_OUT = DERIVED / "recombination_annotated_with_dates.tsv"
SUMMARY_OUT = DERIVED / "recombination_events_summary.tsv"
NODE_DATE_OUT = DATA / "descendants_node_dates.tsv"

FULL_DATE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
MONTH_DATE = re.compile(r"^(\d{4})-(\d{2})$")
YEAR_DATE = re.compile(r"^(\d{4})$")


def normalize_date(raw: str) -> str | None:
    """Return YYYY-MM-DD for full/partial dates, else None."""
    raw = raw.strip()
    m = FULL_DATE.match(raw)
    if m:
        return raw
    m = MONTH_DATE.match(raw)
    if m:
        return f"{m.group(1)}-{m.group(2)}-01"
    m = YEAR_DATE.match(raw)
    if m:
        return f"{m.group(1)}-01-01"
    return None


def main() -> None:
    if not RECOMB_FILE.exists():
        sys.exit(f"missing {RECOMB_FILE} - see data/README.md")
    if not DESC_FILE.exists():
        sys.exit(f"missing {DESC_FILE} (669 MB, not distributed) - see data/README.md")

    # 1. Read recombination events and collect unique node ids.
    nodes = set()
    with open(RECOMB_FILE, encoding="utf-8") as f:
        header = f.readline().rstrip("\n")
        for line in f:
            if not line.strip():
                continue
            nodes.add(line.split("\t", 1)[0])

    print(f"unique recombination events (nodes): {len(nodes)}")

    # 2. Stream descendants.tsv and compute the earliest date per node.
    earliest: dict[str, str | None] = {}
    quality: dict[str, str] = {}
    n_desc: dict[str, int] = {}
    n_dated: dict[str, int] = {}
    n_full_dated: dict[str, int] = {}
    found = 0

    with open(DESC_FILE, encoding="utf-8") as f:
        f.readline()  # header
        for line in f:
            node, rest = line.rstrip("\n").split("\t", 1)
            if node not in nodes:
                continue
            found += 1
            best_full = None
            best_any = None
            n_items = n_full = n_any = 0
            for item in rest.split(","):
                item = item.strip()
                if not item:
                    continue
                n_items += 1
                raw = item.rsplit("|", 1)[-1]
                d_full = normalize_date(raw) if FULL_DATE.match(raw.strip()) else None
                d_any = normalize_date(raw)
                if d_full:
                    n_full += 1
                    if best_full is None or d_full < best_full:
                        best_full = d_full
                if d_any:
                    n_any += 1
                    if best_any is None or d_any < best_any:
                        best_any = d_any
            if best_full is not None:
                earliest[node] = best_full
                quality[node] = "full"
            elif best_any is not None:
                earliest[node] = best_any
                quality[node] = "partial_only"
            else:
                earliest[node] = None
                quality[node] = "missing"
            n_desc[node] = n_items
            n_dated[node] = n_any
            n_full_dated[node] = n_full

    print(f"nodes found in descendants.tsv: {found}")
    missing = [n for n in nodes if earliest.get(n) is None]
    if missing:
        print(f"WARNING: no usable date for {len(missing)} nodes: {missing}")

    # 3. Annotate every recombination row with earliest_date.
    rows = []
    with open(RECOMB_FILE, encoding="utf-8") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            node = line.split("\t", 1)[0]
            rows.append(line.rstrip("\n") + "\t" + (earliest.get(node) or ""))

    with open(ANNOTATED_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(header + "\tearliest_date\n")
        f.write("\n".join(rows) + "\n")
    print(f"wrote {ANNOTATED_OUT}: {len(rows)} rows")

    # 4. One-row-per-event summary (confidence = best recomb_parsimony per node).
    best_parsimony: dict[str, int] = defaultdict(int)
    with open(RECOMB_FILE, encoding="utf-8") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            rp = int(p[-1])
            if rp > best_parsimony[p[0]]:
                best_parsimony[p[0]] = rp

    with open(SUMMARY_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(
            "recomb_node_id\tearliest_date\tn_descendants\tn_dated_descendants\t"
            "n_full_date_descendants\tdate_quality\tbest_recomb_parsimony\n"
        )
        for node in sorted(nodes):
            f.write(
                f"{node}\t{earliest.get(node) or ''}\t{n_desc.get(node, 0)}\t"
                f"{n_dated.get(node, 0)}\t{n_full_dated.get(node, 0)}\t"
                f"{quality.get(node, '')}\t{best_parsimony[node]}\n"
            )
    print(f"wrote {SUMMARY_OUT}: {len(nodes)} events")

    # 5. Compact node -> date table used by the R scripts (small, committable).
    with open(NODE_DATE_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("#node_id\tearliest_date\tn_descendants\tn_dated_descendants\t"
                "n_full_date_descendants\tdate_quality\n")
        for node in sorted(nodes):
            f.write(f"{node}\t{earliest.get(node) or ''}\t{n_desc.get(node, 0)}\t"
                    f"{n_dated.get(node, 0)}\t{n_full_dated.get(node, 0)}\t"
                    f"{quality.get(node, '')}\n")
    print(f"wrote {NODE_DATE_OUT}: {len(nodes)} nodes")

    supported = sum(1 for n in nodes if best_parsimony[n] >= RECOMB_PARSIMONY_MIN)
    print(f"events with recomb_parsimony >= {RECOMB_PARSIMONY_MIN}: {supported}")


if __name__ == "__main__":
    main()
