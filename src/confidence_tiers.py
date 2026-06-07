"""
NBLAST-confidence tier analysis (§4.5).
=======================================

Confidence proxy per triplet = agreement between the curated cross-dataset
match and the NBLAST top-1 match, in {0,1,2}:
    (fafb_match == fafb_nblast_match) + (manc_match == manc_nblast_match)

For each top-k% confidence pool we restrict the edge sets to that pool, take
the giant consensus component, and solve MCIS (best-of-5). A monotonic curve
shows the result is not driven by a pocket of low-confidence matches.

Run:  MCIS_DATA_DIR=/path/to/data python src/confidence_tiers.py
Writes results/confidence_tiers.json.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
from run_analysis import build_solver, resolve_data_dir, _best_of, REPO  # noqa


def main():
    data_dir = resolve_data_dir(None)
    print(f"Data directory: {data_dir}")
    solver = build_solver(data_dir, n_seeds=5)
    solver._load()
    trip = solver._triples.copy()

    def agree(col, nbl):
        if col in trip and nbl in trip:
            return (trip[col].astype(str) == trip[nbl].astype(str)).astype(int)
        return np.zeros(len(trip), int)

    conf = (agree("fafb_match", "fafb_nblast_match")
            + agree("manc_match", "manc_nblast_match")).to_numpy()
    # stable tie-break so the ordering is deterministic
    order = np.lexsort((np.arange(len(trip)), -conf))

    rows = []
    for pct in [10, 20, 30, 50, 75, 100]:
        k = max(1, int(round(len(trip) * pct / 100)))
        keep = set(order[:k].tolist())
        be = {(i, j) for i, j in solver._be if i in keep and j in keep}
        fe = {(i, j) for i, j in solver._fe if i in keep and j in keep}
        me = {(i, j) for i, j in solver._me if i in keep and j in keep}
        CG = nx.DiGraph(); CG.add_edges_from(be & fe & me)
        if not CG.nodes:
            rows.append({"pct": pct, "pool": k, "N": 0}); continue
        giant = sorted(max(nx.weakly_connected_components(CG), key=len))
        loc = {g: l for l, g in enumerate(giant)}
        gbe = {(loc[u], loc[v]) for u, v in be if u in loc and v in loc}
        gfe = {(loc[u], loc[v]) for u, v in fe if u in loc and v in loc}
        gme = {(loc[u], loc[v]) for u, v in me if u in loc and v in loc}
        solver._ng = len(giant)
        n = len(_best_of(solver, gbe, gfe, gme, 5))
        rows.append({"pct": pct, "pool": k, "N": n})
        print(f"  top {pct:>3d}% | pool {k:>4d} | N = {n}")

    (REPO / "results").mkdir(exist_ok=True)
    with open(REPO / "results" / "confidence_tiers.json", "w") as f:
        json.dump({"tiers": rows}, f, indent=2)
    print(f"  Wrote {REPO/'results'/'confidence_tiers.json'}")


if __name__ == "__main__":
    main()
