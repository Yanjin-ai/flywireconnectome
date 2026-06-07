"""
mcis-watch — connectome version-QC command-line tool (Line A).
==============================================================

Given a proofreading edit to one connectome, report which published conserved
(all-3) edges it would gain or lose, named by cell type. The edit is supplied
as a small CSV of ops in that connectome's neuron-id space:

    op,source,target
    edge_delete,720575940623122125,720575940609488942
    edge_insert,720575940623122125,720575940611111111
    node_merge,720575940623122125,720575940622222222

(node_merge/node_split are translated to edge ops internally.)

Usage:
    MCIS_DATA_DIR=/path/to/data mcis-watch --dataset FAFB --edits edits.csv
    # or, after `pip install -e .`:  mcis-watch --dataset FAFB --edits edits.csv

It loads the canonical triplet space once, maps the edited neuron ids to indices,
runs the O(|ΔE|) consensus-impact query against the other two connectomes, and
prints the human-readable QC report (also written next to the edits file).
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir  # noqa
from run_analysis import build_solver  # noqa
from incremental_mcis import (impact_query, qc_report,  # noqa
                              node_ops_to_edge_ops)

ID_COL = {"BANC": "root_626", "FAFB": "fafb_match", "MANC": "manc_match"}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mcis-watch",
                                 description="Connectome version-QC: does an edit "
                                             "touch a conserved circuit?")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--dataset", choices=["BANC", "FAFB", "MANC"], required=True)
    ap.add_argument("--edits", required=True, help="CSV with columns op,source,target")
    args = ap.parse_args(argv)

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng

    # neuron-id -> giant-local index, for the chosen connectome
    col = ID_COL[args.dataset]
    id2loc = {str(solver._triples.iloc[solver._gl[loc]][col]): loc
              for loc in range(ng)}
    ct_map = {loc: str(solver._triples.iloc[solver._gl[loc]].get("cell_type", loc))
              for loc in range(ng)}
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    cur = {"BANC": be, "FAFB": fe, "MANC": me}[args.dataset]

    edf = pd.read_csv(args.edits, dtype=str)
    edge_ops, node_ops, skipped = [], [], 0
    for _, r in edf.iterrows():
        op, s, t = r["op"], str(r["source"]), str(r["target"])
        if s not in id2loc or t not in id2loc:
            skipped += 1
            continue
        i, j = id2loc[s], id2loc[t]
        if op in ("edge_insert", "edge_delete"):
            edge_ops.append((op, i, j))
        elif op == "node_merge":
            node_ops.append(("node_merge", i, j))
        elif op == "node_split":
            node_ops.append(("node_split", i, j, [(i, j)]))
    edge_ops += node_ops_to_edge_ops(node_ops, cur)

    impact = impact_query(be, fe, me, args.dataset, edge_ops)
    report = qc_report(impact, ct_map)
    if skipped:
        report += (f"\n({skipped} edited pair(s) involve neurons outside the "
                   "matched-triplet pool — they cannot affect a conserved edge.)")
    print(report)
    out = os.path.splitext(args.edits)[0] + "_qc_report.txt"
    with open(out, "w") as f:
        f.write(report + "\n")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
