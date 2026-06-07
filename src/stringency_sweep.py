"""
Robustness survey: is N stable under connectomic reconstruction error? (F)
==========================================================================

A strict reviewer asks whether the conserved circuit is an artifact of the
particular edge sets, or robust to the proofreading errors every connectome
carries. We perturb each connectome's edges independently — flipping a fraction
p of edges (remove p·|E| real edges = false negatives, add the same number of
random edges = false positives) — and recompute the best-of-k MCIS. A circuit
that degrades gracefully with p is robust; a cliff would mean fragility.

Run:  MCIS_DATA_DIR=/path/to/data python src/stringency_sweep.py
Outputs: results/stringency_sweep.json, figures/figure13_stringency.png
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir  # noqa
from run_analysis import build_solver, _best_of  # noqa


def perturb(edges, ng, p, rng):
    """Flip a fraction p of a directed edge set: remove p·|E| existing edges and
    add the same number of random non-edges (both false negatives and positives)."""
    E = list(edges)
    if not E:
        return set(E)
    k = int(round(p * len(E)))
    keep = set(E)
    if k > 0:
        rm = rng.choice(len(E), size=min(k, len(E)), replace=False)
        for idx in rm:
            keep.discard(E[idx])
        added = 0
        tries = 0
        while added < k and tries < 20 * k:
            i, j = int(rng.integers(ng)), int(rng.integers(ng))
            tries += 1
            if i != j and (i, j) not in keep:
                keep.add((i, j))
                added += 1
    return keep


def main():
    import argparse
    from pathlib import Path
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--ps", type=float, nargs="+", default=[0.0, 0.05, 0.10, 0.20])
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--k", type=int, default=3)
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)

    rows = []
    for p in args.ps:
        ns = []
        for r in range(args.reps if p > 0 else 1):
            rng = np.random.default_rng(1000 * int(p * 100) + r)
            bep = perturb(be, ng, p, rng)
            fep = perturb(fe, ng, p, rng)
            mep = perturb(me, ng, p, rng)
            ns.append(len(_best_of(solver, bep, fep, mep, args.k)))
        ns = np.array(ns)
        rows.append({"p": p, "mean": float(ns.mean()), "std": float(ns.std()),
                     "min": int(ns.min()), "max": int(ns.max()), "n": ns.tolist()})
        print(f"  p={p:.2f}: N = {ns.mean():.1f} ± {ns.std():.1f} "
              f"[{ns.min()},{ns.max()}]")

    out = {"perturbation": "independent edge flip (FN+FP) per connectome",
           "k_multistart": args.k, "reps": args.reps, "tiers": rows,
           "baseline_N": rows[0]["mean"]}
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "stringency_sweep.json", "w") as f:
        json.dump(out, f, indent=2)

    ps = [r["p"] * 100 for r in rows]
    mean = [r["mean"] for r in rows]
    std = [r["std"] for r in rows]
    fig, ax = plt.subplots(figsize=(7, 5)); fig.patch.set_facecolor("white")
    ax.errorbar(ps, mean, yerr=std, marker="o", lw=2, capsize=4, color="#1f77b4")
    ax.set_xlabel("edge perturbation per connectome (% flipped)")
    ax.set_ylabel("MCIS size N (best-of-%d)" % args.k)
    ax.set_title("Robustness of the conserved circuit to reconstruction error",
                 fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure13_stringency.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print("  Wrote results/stringency_sweep.json, figures/figure13_stringency.png")


if __name__ == "__main__":
    main()
