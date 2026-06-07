"""
Incremental MCIS for connectome version QC (Phase 2 / line A).
==============================================================

Problem. FlyWire-style connectomes change continuously (proofreading edits).
Re-running the full MCIS after every version bump is wasteful. Given the old
MCIS and the edge delta ΔE between versions, can we update the MCIS by touching
only the affected part of the graph?

Formalisation. With neurons bijectively labelled across the three connectomes,
the induced subgraph S is mutually isomorphic iff it contains no "disagreement
edge" — a pair present in some but not all three connectomes. So the MCIS is a
Maximum Independent Set (MIS) on the undirected DISAGREEMENT GRAPH D, where
{i,j} ∈ D iff the pair (i,j) or (j,i) is a disagreement (and any node with a
self-disagreement is excluded outright).

Locality theorem (why incremental works).
  MIS is separable over connected components: an MIS of D is the union of an MIS
  of each connected component of D. An edge delta ΔE only changes D-edges
  incident to the endpoints of the changed pairs. Therefore, after a delta, ONLY
  the connected components of (D_old ∪ D_new) that contain a touched endpoint can
  change; every other component's contribution to S is unchanged. Re-solving just
  those components and keeping the rest fixed yields exactly the same solution a
  full re-solve (with the same component-local greedy) would produce — in time
  proportional to the affected components, not the whole graph.

This module provides: the disagreement-graph construction, a component-local
greedy MIS, the connectome "git diff" (edit ops), the incremental updater, and a
synthetic-delta benchmark (incremental == full re-solve, with speedup). A hook is
provided for real CAVE edit-history deltas.

Run:  MCIS_DATA_DIR=/path/to/data python src/incremental_mcis.py
Outputs: results/incremental_benchmark.json, figures/figure11_incremental.png
"""
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir  # noqa
from run_analysis import build_solver  # noqa


# ─────────────────────────── disagreement graph ───────────────────────────
def build_disagreement(be, fe, me, ng):
    """Return (adj, forced_out): undirected disagreement adjacency over node
    indices, and the set of nodes with a self-disagreement (never in S)."""
    union = be | fe | me
    consensus = be & fe & me
    adj = defaultdict(set)
    forced_out = set()
    for (i, j) in union:
        if (i, j) in consensus:
            continue
        if i == j:
            forced_out.add(i)
        else:
            adj[i].add(j)
            adj[j].add(i)
    return adj, forced_out


def greedy_mis_component(comp, adj, seed=0):
    """Component-local greedy MIS: repeatedly drop the highest-D-degree node
    until the induced subgraph is edgeless; survivors are the independent set."""
    rng = np.random.default_rng(seed)
    active = set(comp)
    # local degree within the component
    deg = {v: len(adj[v] & active) for v in active}
    while True:
        # any remaining edge?
        worst, wd = None, -1
        for v in active:
            if deg[v] > wd or (deg[v] == wd and rng.random() < 0.5):
                worst, wd = v, deg[v]
        if wd <= 0:
            break
        active.discard(worst)
        for u in adj[worst] & active:
            deg[u] -= 1
        del deg[worst]
    return active


def connected_components(nodes, adj):
    seen, comps = set(), []
    for s in nodes:
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


def full_solve(adj, all_nodes, forced_out, seed=0):
    """MIS over the whole disagreement graph = union of per-component MIS.
    Isolated nodes (no disagreement) are automatically in S."""
    candidates = [v for v in all_nodes if v not in forced_out]
    comps = connected_components(candidates, adj)
    S = set()
    for comp in comps:
        S |= greedy_mis_component(comp, adj, seed)
    return S


# ─────────────────────────── connectome "git diff" ─────────────────────────
def diff_edges(old_edges, new_edges):
    """Return edit ops between two directed edge sets of one connectome."""
    ins = new_edges - old_edges
    dele = old_edges - new_edges
    return ([("edge_insert", i, j) for (i, j) in ins]
            + [("edge_delete", i, j) for (i, j) in dele])


# ─────────────────────────── incremental update ────────────────────────────
def _ball(touched, adj, radius, candidates):
    """Radius-r BFS ball (over candidate nodes) around the touched set."""
    frontier = set(v for v in touched if v in candidates)
    ball = set(frontier)
    for _ in range(radius):
        nxt = set()
        for v in frontier:
            nxt |= (adj[v] & candidates) - ball
        ball |= nxt
        frontier = nxt
        if not frontier:
            break
    return ball


def incremental_update(S, adj_old, adj_new, touched, all_nodes, forced_out,
                       radius=None, seed=0):
    """Incrementally update the MIS after an edge delta.

    radius=None  → exact: re-solve every connected component touched by the
                   delta (provably identical to a full re-solve).
    radius=r     → bounded: re-solve only the radius-r ball around the touched
                   nodes, with boundary membership fixed (O(|ball|) work,
                   small controlled approximation).
    Returns (S_new, n_affected_nodes).
    """
    cand = set(v for v in all_nodes if v not in forced_out)
    if radius is None:
        affected = set(v for comp in connected_components(
            [v for v in touched if v in cand], adj_new) for v in comp)
    else:
        affected = _ball(touched, adj_new, radius, cand)
    if not affected:
        return (S - forced_out), 0

    # boundary nodes currently in S that sit just outside the ball forbid their
    # neighbours inside the ball from joining S (independence must hold).
    fixed_in = set(v for v in S if v not in affected and v not in forced_out)
    forbidden = set()
    for v in fixed_in:
        forbidden |= (adj_new[v] & affected)

    resolvable = affected - forbidden
    comps = connected_components(list(resolvable), adj_new)
    S_ball = set()
    for comp in comps:
        comp = [v for v in comp if v in resolvable]
        S_ball |= greedy_mis_component(comp, adj_new, seed)

    S_new = (set(v for v in S if v not in affected) | S_ball) - forced_out
    return S_new, len(affected)


def apply_delta(be, fe, me, dataset, ops):
    """Apply edit ops to one dataset's edge set; return new (be,fe,me)."""
    sets = {"BANC": set(be), "FAFB": set(fe), "MANC": set(me)}
    s = sets[dataset]
    for op, i, j in ops:
        if op == "edge_insert":
            s.add((i, j))
        elif op == "edge_delete":
            s.discard((i, j))
    return sets["BANC"], sets["FAFB"], sets["MANC"]


def touched_nodes(ops):
    t = set()
    for _, i, j in ops:
        t.add(i)
        t.add(j)
    return t


def impact_query(be, fe, me, dataset, ops):
    """The practically useful QC primitive: given an edit to ONE connectome,
    report which CONSENSUS (all-3) edges are gained or lost — in O(|ops|) time,
    independent of graph size, by checking only the edited pairs against the
    other two connectomes. This is what a proofreader actually needs ("does my
    edit touch any published conserved circuit?"), and unlike global-MCIS
    maintenance it is genuinely local because the disagreement graph is dense.
    """
    others = {"BANC": (fe, me), "FAFB": (be, me), "MANC": (be, fe)}[dataset]
    cur = {"BANC": be, "FAFB": fe, "MANC": me}[dataset]
    gained, lost, neurons = [], [], set()
    for op, i, j in ops:
        in_others = (i, j) in others[0] and (i, j) in others[1]
        was_present = (i, j) in cur
        if not in_others:
            continue  # can never be a consensus edge regardless of this edit
        if op == "edge_insert" and not was_present:
            gained.append((i, j)); neurons |= {i, j}
        elif op == "edge_delete" and was_present:
            lost.append((i, j)); neurons |= {i, j}
    return {"consensus_gained": gained, "consensus_lost": lost,
            "neurons_affected": sorted(neurons)}


# ─────────────────────────── benchmark ─────────────────────────────────────
def main():
    import argparse
    from pathlib import Path
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--batches", type=int, default=40)
    ap.add_argument("--delta-size", type=int, default=15,
                    help="edges flipped per version bump")
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    all_nodes = list(range(ng))

    adj, forced = build_disagreement(be, fe, me, ng)
    S = full_solve(adj, all_nodes, forced)
    print(f"  Giant component: {ng} nodes | initial MCIS (MIS) N = {len(S)}")

    # The disagreement graph is essentially one giant component, so EXACT
    # component-decomposition re-solves everything (no speedup). We therefore
    # sweep a k-hop BOUNDED updater and measure the speed/accuracy tradeoff,
    # with radius=None (exact) as the correctness anchor.
    radii = [1, 2, 3, None]
    stats = {str(r): {"speedup": [], "affected": [], "exact": 0,
                      "size_err": []} for r in radii}
    n_comp = len(connected_components(
        [v for v in all_nodes if v not in forced], adj))
    print(f"  Disagreement graph: {n_comp} components "
          f"(largest dominates → exact decomposition gives no locality)")

    rng = np.random.default_rng(0)
    union = list(be | fe | me)
    for b in range(args.batches):
        ds = ["BANC", "FAFB", "MANC"][b % 3]
        ops = []
        for _ in range(args.delta_size):
            if rng.random() < 0.5 and union:
                i, j = union[rng.integers(len(union))]
                ops.append(("edge_delete", i, j))
            else:
                i, j = int(rng.integers(ng)), int(rng.integers(ng))
                if i != j:
                    ops.append(("edge_insert", i, j))
        be2, fe2, me2 = apply_delta(be, fe, me, ds, ops)
        adj2, forced2 = build_disagreement(be2, fe2, me2, ng)
        tch = touched_nodes(ops)

        t0 = time.time(); S_full = full_solve(adj2, all_nodes, forced2)
        t_full = time.time() - t0

        for r in radii:
            t0 = time.time()
            S_inc, n_aff = incremental_update(S, adj, adj2, tch, all_nodes,
                                              forced2, radius=r)
            t_inc = time.time() - t0
            st = stats[str(r)]
            if t_inc > 0:
                st["speedup"].append(t_full / t_inc)
            st["affected"].append(n_aff / ng)
            st["size_err"].append(abs(len(S_inc) - len(S_full)))
            if S_inc == S_full:
                st["exact"] += 1
        be, fe, me, adj, S = be2, fe2, me2, adj2, S_full

    # The O(|ΔE|) consensus-impact primitive (graph-size-independent QC query).
    iq_times = []
    rng2 = np.random.default_rng(1)
    for _ in range(200):
        ops = []
        for _ in range(args.delta_size):
            i, j = int(rng2.integers(ng)), int(rng2.integers(ng))
            if i != j:
                ops.append((["edge_insert", "edge_delete"][int(rng2.integers(2))], i, j))
        t0 = time.time(); impact_query(be, fe, me, "FAFB", ops)
        iq_times.append(time.time() - t0)
    iq_mean_us = float(np.mean(iq_times)) * 1e6
    print(f"  Consensus-impact query (O(|ΔE|), {args.delta_size} edits): "
          f"{iq_mean_us:.1f} µs/query (independent of the {ng}-node graph size)")

    out = {"giant_nodes": ng, "batches": args.batches,
           "delta_size": args.delta_size, "n_disagreement_components": n_comp,
           "consensus_impact_query_us": round(iq_mean_us, 1),
           "by_radius": {}}
    print(f"  {'radius':>7} | {'speedup':>8} | {'%resolved':>9} | "
          f"{'exact':>7} | mean|ΔN|")
    for r in radii:
        st = stats[str(r)]
        sp = float(np.mean(st["speedup"])) if st["speedup"] else 1.0
        af = float(np.mean(st["affected"])) * 100
        ex = st["exact"]
        se = float(np.mean(st["size_err"]))
        out["by_radius"][("exact" if r is None else str(r))] = {
            "mean_speedup": round(sp, 2), "mean_pct_resolved": round(af, 1),
            "exact_match": f"{ex}/{args.batches}", "mean_size_error": round(se, 2)}
        print(f"  {('exact' if r is None else r):>7} | {sp:>7.1f}× | {af:>8.1f}% | "
              f"{ex:>3}/{args.batches} | {se:.2f}")

    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "incremental_benchmark.json", "w") as f:
        json.dump(out, f, indent=2)

    # figure: speedup vs radius and accuracy vs radius
    labels = ["r=1", "r=2", "r=3", "exact"]
    sp = [out["by_radius"][k]["mean_speedup"]
          for k in ("1", "2", "3", "exact")]
    se = [out["by_radius"][k]["mean_size_error"]
          for k in ("1", "2", "3", "exact")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor("white")
    axes[0].bar(labels, sp, color="#1f77b4", alpha=0.85)
    axes[0].set_ylabel("mean speedup vs full re-solve")
    axes[0].set_title("A  Bounded-radius speedup", fontweight="bold")
    for i, v in enumerate(sp):
        axes[0].text(i, v, f"{v:.1f}×", ha="center", va="bottom")
    axes[1].bar(labels, se, color="#d62728", alpha=0.85)
    axes[1].set_ylabel("mean |ΔN| vs exact full re-solve")
    axes[1].set_title("B  Approximation cost (0 = exact)", fontweight="bold")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Incremental MCIS for connectome version QC — "
                 "k-hop bounded update (speed/accuracy tradeoff)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure11_incremental.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print("  Wrote results/incremental_benchmark.json, figures/figure11_incremental.png")


if __name__ == "__main__":
    main()
