"""
Spectral relaxation solver for MCIS (Line B completion).
========================================================

MCIS = Maximum Independent Set on the disagreement graph D (a pair of neurons
cannot both be kept if their connection disagrees across connectomes). Exact MIS
is NP-hard; the greedy + ILP solvers are elsewhere in this repo. Here we add a
SPECTRAL solver: it uses the leading eigenvector of D's adjacency to order
neurons by how "central" they are in the constraint graph, then greedily builds
an independent set from the least-constrained outward. Eigen-decomposition is
O(n·#edges) via Lanczos and is far cheaper than ILP at scale.

We validate the spectral solver against greedy and ILP-optimal on sampled
subgraphs (reusing src/exact_ilp.py), reporting solution quality and runtime.

Run:  MCIS_DATA_DIR=/path/to/data python src/spectral_mcis.py
Output: results/spectral_validation.json, figures/figure12_spectral.png
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir  # noqa
from run_analysis import build_solver  # noqa
from incremental_mcis import build_disagreement, greedy_mis_component  # noqa
from exact_ilp import exact_on_subset, greedy_on_subset  # noqa


def spectral_mis(nodes, adj):
    """Spectral MIS heuristic on the disagreement subgraph induced by `nodes`.
    Order nodes by ascending leading-eigenvector centrality (least-constrained
    first), then greedily add if no chosen neighbour."""
    nodes = sorted(nodes)
    idx = {v: k for k, v in enumerate(nodes)}
    n = len(nodes)
    if n == 0:
        return set()
    rows, cols = [], []
    for v in nodes:
        for u in adj[v]:
            if u in idx:
                rows.append(idx[v]); cols.append(idx[u])
    if not rows:
        return set(nodes)  # no constraints → all independent
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = (A + A.T) / 2.0
    try:
        _, vec = spla.eigsh(A, k=1, which="LA")
        centrality = np.abs(vec[:, 0])
    except Exception:
        centrality = np.asarray(A.sum(1)).ravel()  # fallback: degree
    order = np.argsort(centrality)  # least central first
    chosen, blocked = set(), set()
    for k in order:
        v = nodes[k]
        if v in blocked:
            continue
        chosen.add(v)
        blocked |= adj[v]
    return chosen


def main():
    import argparse
    from pathlib import Path
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--tiers", type=int, nargs="+", default=[30, 50, 80, 120])
    ap.add_argument("--per-tier", type=int, default=12)
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    gbe, gfe, gme = solver._gbe, solver._gfe, solver._gme
    adj, _ = build_disagreement(set(gbe), set(gfe), set(gme), ng)
    rng = np.random.default_rng(7)

    rows = []
    for size in args.tiers:
        g_n, s_n, e_n, t_g, t_s, t_e = [], [], [], [], [], []
        for _ in range(args.per_tier):
            sub = set(rng.choice(ng, size=min(size, ng), replace=False).tolist())
            t0 = time.time(); gn = greedy_on_subset(gbe, gfe, gme, sub); t_g.append(time.time()-t0)
            t0 = time.time(); sn = len(spectral_mis(sub, adj)); t_s.append(time.time()-t0)
            t0 = time.time(); en, _ = exact_on_subset(gbe, gfe, gme, sub); t_e.append(time.time()-t0)
            g_n.append(gn); s_n.append(sn); e_n.append(en)
        g_n, s_n, e_n = map(np.array, (g_n, s_n, e_n))
        rows.append({
            "size": size,
            "greedy_pct_opt": round(100 * float((g_n / np.maximum(e_n, 1)).mean()), 1),
            "spectral_pct_opt": round(100 * float((s_n / np.maximum(e_n, 1)).mean()), 1),
            "greedy_ms": round(1e3 * float(np.mean(t_g)), 2),
            "spectral_ms": round(1e3 * float(np.mean(t_s)), 2),
            "ilp_ms": round(1e3 * float(np.mean(t_e)), 2),
        })
        print(f"  {size:>4d} nodes | greedy {rows[-1]['greedy_pct_opt']}% "
              f"spectral {rows[-1]['spectral_pct_opt']}% of ILP | "
              f"t: greedy {rows[-1]['greedy_ms']}ms spectral {rows[-1]['spectral_ms']}ms "
              f"ILP {rows[-1]['ilp_ms']}ms")

    out = {"tiers": rows,
           "summary": {
               "mean_spectral_pct_opt": round(float(np.mean([r["spectral_pct_opt"] for r in rows])), 1),
               "mean_greedy_pct_opt": round(float(np.mean([r["greedy_pct_opt"] for r in rows])), 1),
           }}
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "spectral_validation.json", "w") as f:
        json.dump(out, f, indent=2)

    sizes = [r["size"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor("white")
    axes[0].plot(sizes, [r["greedy_pct_opt"] for r in rows], "o-", label="greedy")
    axes[0].plot(sizes, [r["spectral_pct_opt"] for r in rows], "s-", label="spectral")
    axes[0].axhline(100, color="#999", ls="--", label="ILP optimum")
    axes[0].set_xlabel("subgraph size"); axes[0].set_ylabel("% of ILP optimum")
    axes[0].legend(); axes[0].set_title("A  Solution quality vs exact ILP", fontweight="bold")
    axes[1].plot(sizes, [r["greedy_ms"] for r in rows], "o-", label="greedy")
    axes[1].plot(sizes, [r["spectral_ms"] for r in rows], "s-", label="spectral")
    axes[1].plot(sizes, [r["ilp_ms"] for r in rows], "^-", label="ILP")
    axes[1].set_yscale("log"); axes[1].set_xlabel("subgraph size")
    axes[1].set_ylabel("solve time (ms, log)"); axes[1].legend()
    axes[1].set_title("B  Runtime", fontweight="bold")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Spectral MCIS relaxation — quality/speed vs greedy and ILP",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure12_spectral.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print(f"  Mean spectral quality: {out['summary']['mean_spectral_pct_opt']}% of ILP "
          f"(greedy {out['summary']['mean_greedy_pct_opt']}%)")
    print("  Wrote results/spectral_validation.json, figures/figure12_spectral.png")


if __name__ == "__main__":
    main()
