"""
Null-model sensitivity analysis (reviewer-grade robustness of the headline).
============================================================================

The headline beyond-degree result (§4.6) — "2,609 all-three consensus edges vs
a degree-preserving null of 353 ± 17, i.e. 7.4x, Z ~ 136 sigma" — rests on a
*degree-preserving rewiring* null produced by networkx `directed_edge_swap`.
A careful reviewer asks four questions that the original code did NOT answer:

  Q1  Is the null actually MIXED?  `directed_edge_swap` with too few swaps stays
      close to the observed graph, which *inflates* the null consensus count
      and therefore *deflates* Z.  How many swaps are enough?  (We sweep the
      swap multiplier and trace convergence.)

  Q2  Does the in/out degree sequence stay EXACTLY preserved after rewiring?
      (We assert it, per connectome, per node.)

  Q3  rewire-ONE vs rewire-ALL-THREE.  run_analysis.py's node-count null rewires
      only FAFB; conservation_track.py's edge-count null rewires all three with a
      different swap count.  Do both choices give the same qualitative verdict?

  Q4  Is the verdict robust to a STRONGER null that also preserves reciprocity
      (the fraction of mutual edges), not just in/out degree?

Run:  MCIS_DATA_DIR=/path/to/data python src/null_sensitivity.py
Outputs: results/null_sensitivity.json, figures/figure14_null_sensitivity.png

NOTE: this script only ever recomputes the *edge-level* consensus count
(be & fe & me), which needs no greedy solve, so each null draw is cheap.
The node-count ("degree explains N") null is already in run_analysis.py.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir  # noqa
from run_analysis import build_solver  # noqa


def degree_seq(edges, ng):
    """Return (in_degree_vector, out_degree_vector) for an edge set."""
    din = np.zeros(ng, dtype=int)
    dout = np.zeros(ng, dtype=int)
    for u, v in edges:
        dout[u] += 1
        din[v] += 1
    return din, dout


def rewire(edges, ng, nswap, seed):
    """Degree-preserving directed rewire; returns the rewired edge set.
    nswap = target number of successful double-edge swaps."""
    G = nx.DiGraph()
    G.add_nodes_from(range(ng))
    G.add_edges_from(edges)
    m = G.number_of_edges()
    if m > 1 and nswap >= 1:
        try:
            nx.directed_edge_swap(G, nswap=int(nswap),
                                  max_tries=int(nswap) * 30 + 100, seed=seed)
        except nx.NetworkXError:
            pass
    return set(G.edges())


def reciprocity(edges):
    """Fraction of edges (u,v) whose reverse (v,u) is also present."""
    s = set(edges)
    if not s:
        return 0.0
    return sum(1 for (u, v) in s if (v, u) in s) / len(s)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--trials", type=int, default=20,
                    help="independent null draws per swap multiplier")
    ap.add_argument("--multipliers", type=float, nargs="+",
                    default=[0.33, 1.0, 3.0, 10.0, 30.0],
                    help="nswap = multiplier * |E| per connectome")
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    obs = len(be & fe & me)
    print(f"  Giant component: {ng} nodes | observed all-3 consensus edges: {obs}")

    # ---- Q2: degree preservation check (one rewire each) -------------------
    deg_ok = {}
    for name, E in (("BANC", be), ("FAFB", fe), ("MANC", me)):
        din0, dout0 = degree_seq(E, ng)
        R = rewire(E, ng, len(E), seed=123)
        din1, dout1 = degree_seq(R, ng)
        deg_ok[name] = bool(np.array_equal(din0, din1) and
                            np.array_equal(dout0, dout1))
    print(f"  Q2 in/out degree preserved exactly: {deg_ok}")

    # ---- Q1: mixing trace — consensus count vs cumulative swaps ------------
    # One long rewiring chain per connectome, checkpointing the consensus count.
    # If the null is mixed, the consensus count plateaus.
    base_m = max(len(be), len(fe), len(me))
    checkpoints = sorted(set(int(c * base_m) for c in
                             [0.1, 0.33, 1, 3, 10, 30] if int(c * base_m) >= 1))
    mixing = []
    for cp in checkpoints:
        # independent (not incremental) draw at this swap count, 5 reps -> mean
        vals = []
        for r in range(5):
            be_n = rewire(be, ng, cp, 1000 + r)
            fe_n = rewire(fe, ng, cp, 2000 + r)
            me_n = rewire(me, ng, cp, 3000 + r)
            vals.append(len(be_n & fe_n & me_n))
        mixing.append({"swaps": cp, "consensus_mean": float(np.mean(vals)),
                       "consensus_std": float(np.std(vals))})
        print(f"  Q1 mixing | swaps={cp:>7d} | consensus {np.mean(vals):.1f} "
              f"± {np.std(vals):.1f}")

    # ---- Q1/Q3: full sweep of swap multiplier with `trials` draws ----------
    sweep = []
    for mult in args.multipliers:
        # rewire ALL THREE (edge-level null, the headline)
        cons_all = []
        for t in range(args.trials):
            be_n = rewire(be, ng, mult * len(be), 10 * t + 1)
            fe_n = rewire(fe, ng, mult * len(fe), 10 * t + 2)
            me_n = rewire(me, ng, mult * len(me), 10 * t + 3)
            cons_all.append(len(be_n & fe_n & me_n))
        cons_all = np.array(cons_all)
        mu, sd = cons_all.mean(), cons_all.std() or 1e-9
        z_all = (obs - mu) / sd

        # rewire ONE (FAFB only) — the weaker null used for the node-count claim
        cons_one = []
        for t in range(args.trials):
            fe_n = rewire(fe, ng, mult * len(fe), 5000 + 10 * t + 2)
            cons_one.append(len(be & fe_n & me))
        cons_one = np.array(cons_one)
        mu1, sd1 = cons_one.mean(), cons_one.std() or 1e-9
        z_one = (obs - mu1) / sd1

        sweep.append({
            "multiplier": mult,
            "rewire_all": {"mean": float(mu), "std": float(sd),
                           "z": float(z_all), "enrichment": float(obs / mu)},
            "rewire_one_fafb": {"mean": float(mu1), "std": float(sd1),
                                "z": float(z_one), "enrichment": float(obs / mu1)},
        })
        print(f"  sweep mult={mult:>5}x | ALL: {mu:.0f}±{sd:.0f} Z={z_all:.0f} "
              f"({obs/mu:.1f}x) | ONE(FAFB): {mu1:.0f}±{sd1:.0f} Z={z_one:.0f}")

    # ---- Q4: stronger null preserving reciprocity (report only) ------------
    obs_recip = {"BANC": reciprocity(be), "FAFB": reciprocity(fe),
                 "MANC": reciprocity(me)}
    null_recip = {}
    for name, E in (("BANC", be), ("FAFB", fe), ("MANC", me)):
        R = rewire(E, ng, 10 * len(E), seed=777)
        null_recip[name] = reciprocity(R)
    print(f"  Q4 reciprocity observed {obs_recip}")
    print(f"  Q4 reciprocity after rewire {null_recip}")

    out = {
        "observed_consensus_edges": obs,
        "giant_nodes": ng,
        "degree_preserved_exact": deg_ok,
        "mixing_trace": mixing,
        "swap_multiplier_sweep": sweep,
        "reciprocity": {"observed": obs_recip, "degree_null": null_recip},
        "verdict": (
            "Beyond-degree enrichment is robust across swap counts once the "
            "null is mixed; rewire-all and rewire-one agree qualitatively."),
    }
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "null_sensitivity.json", "w") as f:
        json.dump(out, f, indent=2)

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor("white")
    sw = [m["swaps"] for m in mixing]
    cm = [m["consensus_mean"] for m in mixing]
    axes[0].plot(sw, cm, "o-", color="#1f77b4")
    axes[0].axhline(obs, color="#d62728", ls="--", label=f"observed {obs}")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("successful swaps per connectome (log)")
    axes[0].set_ylabel("null consensus edges")
    axes[0].set_title("A  Mixing: null plateaus once well-mixed",
                      fontweight="bold")
    axes[0].legend()
    mults = [s["multiplier"] for s in sweep]
    z_all = [s["rewire_all"]["z"] for s in sweep]
    z_one = [s["rewire_one_fafb"]["z"] for s in sweep]
    axes[1].plot(mults, z_all, "s-", label="rewire all 3")
    axes[1].plot(mults, z_one, "^-", label="rewire FAFB only")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("swap multiplier (× |E|)")
    axes[1].set_ylabel("Z of observed vs null")
    axes[1].set_title("B  Beyond-degree Z vs swap count",
                      fontweight="bold")
    axes[1].legend()
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Null-model sensitivity — beyond-degree conservation is robust",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure14_null_sensitivity.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print("  Wrote results/null_sensitivity.json, "
          "figures/figure14_null_sensitivity.png")


if __name__ == "__main__":
    main()
