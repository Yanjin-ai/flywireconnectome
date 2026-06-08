"""
Canonical analysis driver — single source of truth for all reported numbers.
============================================================================

Runs the MCIS pipeline (MCISSolver: 987-node giant consensus component,
GMIN + (1,2)-swap multi-start = Maximum Independent Set on the disagreement
graph) and all null models from ONE place, so that science.md / README quote
exactly what this script produces.

Data directory resolution (in order):
    1. --data-dir CLI argument
    2. MCIS_DATA_DIR environment variable
    3. ./data/   (relative to repo root)

Expected files in the data directory:
    banc_meta.feather
    banc_626_edge_list.csv   (or 'banc_626_edge_list (2).csv')
    fafb_783_edge_list.csv
    manc_1.2.1_edge_list.csv

Outputs (written into the repo, not the data dir):
    network.csv
    network_enriched.csv
    results/canonical_results.json
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
from mcis_connectome.solver import MCISSolver, MCISResult  # noqa: E402
from mcis_paths import repo_root  # noqa: E402

REPO = Path(repo_root())


def resolve_data_dir(cli_arg):
    for cand in (cli_arg, os.environ.get("MCIS_DATA_DIR"), str(REPO / "data")):
        if cand and Path(cand).is_dir():
            return Path(cand)
    raise SystemExit(
        "Could not locate data directory. Pass --data-dir, set MCIS_DATA_DIR, "
        "or place files under ./data/."
    )


def find_edge_list(data_dir: Path, *names):
    for n in names:
        p = data_dir / n
        if p.exists():
            return str(p)
    raise SystemExit(f"None of {names} found in {data_dir}")


def build_solver(data_dir: Path, n_seeds: int) -> MCISSolver:
    return MCISSolver(
        edge_lists={
            "BANC": find_edge_list(data_dir, "banc_626_edge_list.csv",
                                   "banc_626_edge_list (2).csv"),
            "FAFB": find_edge_list(data_dir, "fafb_783_edge_list.csv"),
            "MANC": find_edge_list(data_dir, "manc_1.2.1_edge_list.csv"),
        },
        triplets_path=str(data_dir / "banc_meta.feather"),
        n_seeds=n_seeds,
    )


def _best_of(solver, gbe, gfe, gme, k, seed0=0):
    """Best MCIS over k restarts of the production GMIN + (1,2)-swap solver."""
    best_set, _ = solver._best_mcis(gbe, gfe, gme, solver._ng,
                                    restarts=k, seed0=seed0)
    return best_set


def seed_distribution(solver: MCISSolver, n_seeds: int):
    """Run the production GMIN + (1,2)-swap solver for n_seeds restarts; return
    the per-restart size distribution and the best node-set found."""
    best_set, sizes = solver._best_mcis(solver._gbe, solver._gfe, solver._gme,
                                        solver._ng, restarts=n_seeds)
    return np.array(sizes), best_set


def correspondence_shuffle_null(solver, n_trials, k=200):
    """Permute the FAFB *and* MANC triplet mapping (full correspondence
    destruction). Each trial = best-of-k multi-start, matching the real
    procedure. Returns null sizes array."""
    ng = solver._ng
    gfe, gme, null = list(solver._gfe), list(solver._gme), []
    for t in range(n_trials):
        rng = np.random.default_rng(1000 + t)
        pf, pm = rng.permutation(ng), rng.permutation(ng)
        gfe_n = {(pf[i], pf[j]) for i, j in gfe}
        gme_n = {(pm[i], pm[j]) for i, j in gme}
        null.append(len(_best_of(solver, solver._gbe, gfe_n, gme_n, k)))
    return np.array(null)


def degree_preserving_null(solver, n_trials, k=200, swap_mult=10):
    """Rewire the FAFB edge set preserving in/out degree per node, then rerun
    (best-of-k). Tests whether degree sequence alone explains achievable N.

    swap_mult * |E| double-edge swaps. The original m//3 (~0.33x) was
    UNDER-MIXED (null 101.0 ~ real, overstating "degree explains all of N");
    well-mixed (>=3x|E|, see src/null_sensitivity.py) the null drops to ~95.5,
    so degree explains MOST (~95%) but not all of N — neuron identity adds the
    remaining ~5 neurons. Default is now a well-mixed 10x|E|."""
    ng = solver._ng
    null = []
    for t in range(n_trials):
        G = nx.DiGraph()
        G.add_nodes_from(range(ng))
        G.add_edges_from(solver._gfe)
        n_e = G.number_of_edges()
        try:
            nx.directed_edge_swap(G, nswap=max(1, swap_mult * n_e),
                                  max_tries=swap_mult * n_e * 30 + 100,
                                  seed=2000 + t)
        except nx.NetworkXError:
            pass  # too few edges to swap; use as-is
        gfe_n = set(G.edges())
        null.append(len(_best_of(solver, solver._gbe, gfe_n, solver._gme, k)))
    return np.array(null)


def centrality_test(solver, best_set, n_perm=1000):
    """Betweenness of circuit vs non-circuit nodes on the consensus graph,
    with a label-permutation p-value."""
    CG = nx.DiGraph()
    CG.add_edges_from(solver._gbe & solver._gfe & solver._gme)
    CG.add_nodes_from(range(solver._ng))
    bc = nx.betweenness_centrality(CG)
    vals = np.array([bc[i] for i in range(solver._ng)])
    in_circ = np.array([i in best_set for i in range(solver._ng)])
    obs = vals[in_circ].mean() - vals[~in_circ].mean()
    rng = np.random.default_rng(42)
    n_c = int(in_circ.sum())
    count = 0
    for _ in range(n_perm):
        idx = rng.permutation(solver._ng)[:n_c]
        mask = np.zeros(solver._ng, bool)
        mask[idx] = True
        if (vals[mask].mean() - vals[~mask].mean()) <= obs:
            count += 1
    return {
        "circuit_mean": float(vals[in_circ].mean()),
        "noncircuit_mean": float(vals[~in_circ].mean()),
        "observed_diff": float(obs),
        "p_value": (count + 1) / (n_perm + 1),
        "n_perm": n_perm,
    }


ENRICHED_COLS = [
    "root_626", "fafb_match", "manc_match", "super_class", "cell_type",
    "sexually_dimorphic", "hemilineage", "neurotransmitter_predicted",
    "cns_network", "fafb_cell_type", "manc_cell_type", "side",
]


def write_enriched(triples_df, path):
    """Write the curated enriched circuit table (stable schema consumed by
    visualize.py and derived_stats.py)."""
    df = triples_df.copy()
    out = pd.DataFrame({
        "BANC": df["root_626"].astype(str),
        "FAFB": df["fafb_match"].astype(str),
        "MANC": df["manc_match"].astype(str),
    })
    for c in ENRICHED_COLS:
        out[c] = df[c] if c in df.columns else pd.NA
    out.to_csv(path, index=False)


def main():
    ap = argparse.ArgumentParser(description="Canonical MCIS analysis driver.")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--seeds", type=int, default=3000,
                    help="GMIN+2-swap restarts for the production circuit")
    ap.add_argument("--null-trials", type=int, default=30)
    args = ap.parse_args()

    data_dir = resolve_data_dir(args.data_dir)
    print(f"Data directory: {data_dir}")
    (REPO / "results").mkdir(exist_ok=True)

    solver = build_solver(data_dir, args.seeds)
    solver._load()

    t0 = time.time()
    sizes, best_set = seed_distribution(solver, args.seeds)
    print(f"\n  {args.seeds}-seed distribution: "
          f"N = {sizes.mean():.1f} ± {sizes.std():.1f}, "
          f"range [{sizes.min()}, {sizes.max()}]  ({time.time()-t0:.1f}s)")
    best_n = int(sizes.max())

    # All nulls use the SAME procedure as the real search (best-of-k multi-start)
    # so comparisons are apples-to-apples. Z is reported against the real best N.
    def z(arr):
        return (best_n - arr.mean()) / arr.std()

    cs = correspondence_shuffle_null(solver, args.null_trials, k=200)
    print(f"  Correspondence-shuffle null ({args.null_trials}×best-of-200): "
          f"{cs.mean():.1f} ± {cs.std():.1f}; Z = {z(cs):.1f}σ")

    dp = degree_preserving_null(solver, args.null_trials, k=200)
    print(f"  Degree-preserving null ({args.null_trials}×best-of-200): "
          f"{dp.mean():.1f} ± {dp.std():.1f}; Z = {z(dp):.1f}σ")

    consensus_edges = len(solver._gbe & solver._gfe & solver._gme)

    # Write canonical network.csv from the best node-set
    global_idx = [solver._gl[i] for i in sorted(best_set)]
    edge_sets = [
        {(i, j) for i, j in solver._gbe if i in best_set and j in best_set},
        {(i, j) for i, j in solver._gfe if i in best_set and j in best_set},
        {(i, j) for i, j in solver._gme if i in best_set and j in best_set},
    ]
    result = MCISResult(best_set,
                        solver._triples.iloc[global_idx].reset_index(drop=True),
                        edge_sets, list(solver.edge_lists.keys()))
    result.to_csv(str(REPO / "network.csv"))
    write_enriched(result.triples, REPO / "network_enriched.csv")
    print("\n" + result.summary())

    # conserved directed edges as cell-type pairs (for the viz / explorer)
    ct = {i: solver._triples.iloc[solver._gl[i]].get("cell_type", "?")
          for i in best_set}
    conserved_edge_celltypes = sorted(
        {(str(ct[i]), str(ct[j])) for (i, j) in edge_sets[0]})

    cent = centrality_test(solver, best_set)
    print(f"  Centrality: circuit betweenness {cent['circuit_mean']:.2e} vs "
          f"non-circuit {cent['noncircuit_mean']:.2e} (p = {cent['p_value']:.3f})")

    canonical = {
        "N_reported": best_n,
        "n_conserved_edges": result.n_edges,
        "isomorphic": bool(result.is_isomorphic),
        "conserved_edge_celltypes": [list(e) for e in conserved_edge_celltypes],
        "consensus_edges_giant": consensus_edges,
        "seed_distribution": {
            "n_seeds": int(args.seeds),
            "mean": float(sizes.mean()), "std": float(sizes.std()),
            "min": int(sizes.min()), "max": int(sizes.max()),
            "sizes": sizes.tolist(),
        },
        "correspondence_shuffle_null": {
            "n_trials": int(args.null_trials), "per_trial": "best-of-200",
            "mean": float(cs.mean()), "std": float(cs.std()),
            "z_vs_best": float(z(cs)), "sizes": cs.tolist(),
        },
        "degree_preserving_null": {
            "n_trials": int(args.null_trials), "per_trial": "best-of-200",
            "mean": float(dp.mean()), "std": float(dp.std()),
            "z_vs_best": float(z(dp)), "sizes": dp.tolist(),
        },
        "centrality": cent,
        "n_bounds": {
            "matched_triplets": int(len(solver._meta)),
            "in_all3_edge_lists": int(len(solver._triples)),
            "giant_component": int(solver._ng),
            "consensus_edges": consensus_edges,
            "best_N": best_n,
        },
    }
    out = REPO / "results" / "canonical_results.json"
    with open(out, "w") as f:
        json.dump(canonical, f, indent=2)
    print(f"\n  Wrote {out}")


if __name__ == "__main__":
    main()
