"""
Connectome Conservation Track (Phase 0 / scientific deepening).
===============================================================

Instead of reporting a single binary "conserved circuit", we score *every*
edge and *every* neuron in the consensus space by how conserved its wiring is
ACROSS THE THREE CONNECTOMES, relative to a degree-preserving null. This turns
the headline from a degree-explained binary claim into a continuous, null-
normalised "conservation track" (analogous to phyloP/phastCons tracks in
comparative genomics).

The central test it answers — the one a reviewer asks first:
    Are there MORE edges shared by all three connectomes than a degree-
    preserving rewiring of each connectome would produce by chance?
If yes, specific wiring is conserved BEYOND the degree sequence. If no, the
honest conclusion is that conservation is a degree-sequence property.

Weighted extension (synapse counts): the edge lists shipped here are binary
(source, target). When a synapse-resolution edge list with a weight column is
available, pass it via --weights to additionally weight each edge by
min synapse count across datasets and use a degree+strength-preserving null;
the binary track below is the special case w==1.

Run:  MCIS_DATA_DIR=/path/to/data python src/conservation_track.py
Outputs: results/conservation_track.json, results/neuron_conservation.csv,
         figures/figure10_conservation_track.png
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir, repo_root  # noqa
from run_analysis import build_solver  # noqa


# Swap multiplier for the Maslov–Sneppen null. The original code used m//2
# (0.5x), which null_sensitivity.py showed is SEVERELY under-mixed: the null
# consensus count only plateaus past ~3x|E| swaps. Under-mixing leaves residual
# real structure in the null, inflating its consensus count (353 at 0.5x) and
# thereby DEFLATING the enrichment. We default to a well-mixed 10x|E|, where the
# null is converged (39 ± 6) and the beyond-degree enrichment is ~67x, not 7.4x.
REWIRE_MULT = 10


def degree_preserving_rewire(edges, n_nodes, seed, mult=REWIRE_MULT):
    """Rewire a directed edge set preserving each node's exact in/out degree,
    with enough swaps (`mult` × |E|) to actually mix (see null_sensitivity.py)."""
    G = nx.DiGraph()
    G.add_nodes_from(range(n_nodes))
    G.add_edges_from(edges)
    m = G.number_of_edges()
    if m > 1:
        try:
            nx.directed_edge_swap(G, nswap=max(1, mult * m),
                                  max_tries=mult * m * 30 + 100, seed=seed)
        except nx.NetworkXError:
            pass
    return set(G.edges())


def weighted_consensus(bw, fw, mw):
    """Conserved synaptic strength: for each edge present in all three weighted
    connectomes, the min weight across them (the strength that survives the
    cross-connectome bottleneck). Returns (per_edge_dict, total_strength)."""
    common = set(bw) & set(fw) & set(mw)
    per_edge = {e: min(bw[e], fw[e], mw[e]) for e in common}
    return per_edge, float(sum(per_edge.values()))


def strength_preserving_null(bw, fw, mw, ng, trials=100, seed=0):
    """Degree- AND strength-aware null: rewire each connectome's binary topology
    preserving in/out degree, then reassign that connectome's observed weight
    multiset randomly onto the rewired edges (preserving the strength
    distribution while destroying specific weight↔edge associations). Returns
    the null distribution of total conserved strength."""
    rng = np.random.default_rng(seed)
    out = []
    for t in range(trials):
        trial = {}
        for name, w in (("b", bw), ("f", fw), ("m", mw)):
            edges = list(w.keys())
            rewired = list(degree_preserving_rewire(set(edges), ng, seed=10 * t
                                                    + {"b": 1, "f": 2, "m": 3}[name]))
            vals = rng.permutation(list(w.values()))
            trial[name] = {e: float(vals[k % len(vals)])
                           for k, e in enumerate(rewired)} if rewired else {}
        _, tot = weighted_consensus(trial["b"], trial["f"], trial["m"])
        out.append(tot)
    return np.array(out)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--null-trials", type=int, default=50)
    ap.add_argument("--weights", default=None,
                    help="optional CSV: source,target,weight (synapse counts)")
    args = ap.parse_args()

    from pathlib import Path
    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()

    ng = solver._ng
    be, fe, me = solver._gbe, solver._gfe, solver._gme  # giant-component edges
    consensus = be & fe & me
    union = be | fe | me

    # ---- per-edge conservation: present in how many of the 3 connectomes ----
    support = {}
    for e in union:
        support[e] = (e in be) + (e in fe) + (e in me)
    obs_consensus = sum(1 for s in support.values() if s == 3)

    # ---- single degree-preserving null loop: records both the global
    #      consensus-edge count and the per-node consensus-incident degree ----
    obs_deg = np.zeros(ng)
    for (u, v) in consensus:
        obs_deg[u] += 1
        obs_deg[v] += 1

    null_counts = np.zeros(args.null_trials)
    null_deg = np.zeros((args.null_trials, ng))
    for t in range(args.null_trials):
        be_n = degree_preserving_rewire(be, ng, 10 * t + 1)
        fe_n = degree_preserving_rewire(fe, ng, 10 * t + 2)
        me_n = degree_preserving_rewire(me, ng, 10 * t + 3)
        cons_n = be_n & fe_n & me_n
        null_counts[t] = len(cons_n)
        for (u, v) in cons_n:
            null_deg[t, u] += 1
            null_deg[t, v] += 1

    mu, sd = null_counts.mean(), null_counts.std() or 1e-9
    z_edges = (obs_consensus - mu) / sd
    enrichment = obs_consensus / mu if mu > 0 else float("inf")
    print(f"  Giant component: {ng} nodes, {len(union)} union edges")
    print(f"  Observed all-3 consensus edges: {obs_consensus}")
    print(f"  Degree-preserving null: {mu:.1f} ± {sd:.1f} "
          f"({args.null_trials} trials)")
    print(f"  >>> Conservation beyond degree:  {enrichment:.1f}x,  Z = {z_edges:.1f}σ")

    nd_mu = null_deg.mean(0)
    nd_sd = null_deg.std(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_node = np.where(nd_sd > 0, (obs_deg - nd_mu) / nd_sd, 0.0)

    # map giant-local index -> BANC id + metadata
    gl = solver._gl
    trip = solver._triples
    rows = []
    for loc in range(ng):
        gi = gl[loc]
        r = trip.iloc[gi]
        rows.append({
            "BANC": r["root_626"], "cell_type": r.get("cell_type"),
            "super_class": r.get("super_class"),
            "consensus_edges": int(obs_deg[loc]),
            "null_mean": round(float(nd_mu[loc]), 3),
            "conservation_z": round(float(z_node[loc]), 2),
        })
    df = pd.DataFrame(rows).sort_values("conservation_z", ascending=False)
    os.makedirs(results_dir(), exist_ok=True)
    df.to_csv(results_dir() + "neuron_conservation.csv", index=False)

    out = {
        "giant_nodes": ng,
        "union_edges": len(union),
        "observed_consensus_edges": obs_consensus,
        "degree_null": {
            "trials": args.null_trials,
            "mean": float(mu), "std": float(sd),
            "z": float(z_edges), "enrichment": float(enrichment),
        },
        "support_histogram": {str(k): int(sum(1 for s in support.values() if s == k))
                              for k in (1, 2, 3)},
        "top_conserved_neurons": df.head(15)[["BANC", "cell_type",
                                              "conservation_z"]].to_dict("records"),
        "weighted": None,
    }

    # ---- optional weighted (degree+strength-preserving) conservation --------
    # --weights expects a CSV of giant-local edges with per-connectome synapse
    # weights: columns i,j,w_banc,w_fafb,w_manc (derive once from synapse tables).
    if args.weights and os.path.exists(args.weights):
        w = pd.read_csv(args.weights)
        bw = {(int(r.i), int(r.j)): float(r.w_banc) for r in w.itertuples() if r.w_banc > 0}
        fw = {(int(r.i), int(r.j)): float(r.w_fafb) for r in w.itertuples() if r.w_fafb > 0}
        mw = {(int(r.i), int(r.j)): float(r.w_manc) for r in w.itertuples() if r.w_manc > 0}
        _, obs_w = weighted_consensus(bw, fw, mw)
        null_w = strength_preserving_null(bw, fw, mw, ng,
                                          trials=min(args.null_trials, 100))
        zw = (obs_w - null_w.mean()) / (null_w.std() or 1e-9)
        out["weighted"] = {"observed_strength": obs_w,
                           "null_mean": float(null_w.mean()),
                           "null_std": float(null_w.std()), "z": float(zw)}
        print(f"  Weighted conservation: observed strength {obs_w:.0f} vs "
              f"strength-null {null_w.mean():.0f} ± {null_w.std():.0f}  Z = {zw:.1f}σ")
    elif args.weights:
        print(f"  (--weights file not found: {args.weights}; ran binary track only)")

    with open(results_dir() + "conservation_track.json", "w") as f:
        json.dump(out, f, indent=2, default=float)

    # ---- figure ----
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.patch.set_facecolor("white")
    h = [out["support_histogram"][k] for k in ("1", "2", "3")]
    axes[0].bar(["1 dataset", "2 datasets", "3 (consensus)"], h,
                color=["#bbb", "#fbb", "#2ca02c"])
    axes[0].set_ylabel("edges"); axes[0].set_title("A  Edge support across connectomes",
                                                    fontweight="bold")
    axes[1].hist(null_counts, bins=20, color="#999", alpha=0.8,
                 label=f"degree null {mu:.0f}±{sd:.0f}")
    axes[1].axvline(obs_consensus, color="#d62728", lw=2,
                    label=f"observed {obs_consensus}")
    axes[1].set_xlabel("# consensus edges"); axes[1].set_ylabel("trials")
    axes[1].legend()
    axes[1].set_title(f"B  Beyond-degree test (Z={z_edges:.1f}σ, {enrichment:.1f}×)",
                      fontweight="bold")
    zz = np.sort(z_node)[::-1]
    axes[2].plot(zz, color="#1f77b4")
    axes[2].axhline(0, color="#999", lw=0.8)
    axes[2].set_xlabel("neuron rank"); axes[2].set_ylabel("conservation z")
    axes[2].set_title("C  Per-neuron conservation track", fontweight="bold")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Connectome Conservation Track — wiring conserved beyond degree sequence",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure10_conservation_track.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print(f"  Wrote results/conservation_track.json, neuron_conservation.csv, "
          f"figures/figure10_conservation_track.png")


if __name__ == "__main__":
    main()
