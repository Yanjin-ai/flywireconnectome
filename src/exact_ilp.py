"""
Exact MCIS validation via Integer Linear Programming (PuLP / CBC).
==================================================================

Provides a certified optimality benchmark for the greedy + exhaustive
expansion heuristic, on randomly sampled induced subgraphs of the 987-node
consensus component.

Formulation
-----------
With nodes bijectively labelled across the three connectomes, the induced
subgraph S is mutually isomorphic iff it contains NO "disagreement edge"
— an ordered pair (i, j) present in some but not all three datasets. So:

    maximise   sum_v x_v
    subject to x_i + x_j <= 1   for every disagreement edge (i, j)
               x_v in {0, 1}

This is Maximum Independent Set on the (undirected) disagreement graph,
solved to provable optimality with CBC for small instances.

Usage:
    python src/exact_ilp.py --data-dir /path/to/data \
        --tiers 20 30 40 50 --per-tier 15 15 10 10
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pulp

sys.path.insert(0, os.path.dirname(__file__))
from run_analysis import build_solver, resolve_data_dir, REPO  # noqa: E402
from incremental_mcis import build_disagreement  # noqa: E402


def make_sampler(adj, ng, kind):
    """Return a function rng,size -> set of local node indices.

    The greedy optimality gap depends critically on WHERE we sample. A node is
    "hard" iff it sits in a dense disagreement neighbourhood — that is exactly
    where max-degree greedy peels too aggressively (science.md §3.4). Uniform
    sampling draws mostly sparse, easy subgraphs and therefore reports an
    OPTIMISTIC gap. We expose three schemes:

      uniform            — original: each node equally likely (easy regime)
      degree_stratified  — P(node) ∝ 1 + disagreement-degree (over-weights the
                           constrained core, a harder regime)
      disagreement_ego   — grow a connected ball in the disagreement graph from
                           a high-degree seed (the worst case for greedy)
    """
    deg = np.array([len(adj.get(v, ())) for v in range(ng)], dtype=float)

    def uniform(rng, size):
        return set(rng.choice(ng, size=min(size, ng), replace=False).tolist())

    def degree_stratified(rng, size):
        w = 1.0 + deg
        p = w / w.sum()
        return set(rng.choice(ng, size=min(size, ng), replace=False,
                              p=p).tolist())

    def disagreement_ego(rng, size):
        # BFS a connected ball in D from a seed chosen ∝ disagreement degree
        if deg.sum() == 0:
            return uniform(rng, size)
        seed = int(rng.choice(ng, p=(deg / deg.sum())))
        chosen, frontier = {seed}, [seed]
        while frontier and len(chosen) < size:
            nxt = []
            rng.shuffle(frontier)
            for u in frontier:
                for w_ in adj.get(u, ()):
                    if w_ not in chosen:
                        chosen.add(w_)
                        nxt.append(w_)
                        if len(chosen) >= size:
                            break
                if len(chosen) >= size:
                    break
            frontier = nxt
        # top up with random nodes if the ball was too small
        if len(chosen) < size:
            extra = [v for v in rng.permutation(ng).tolist() if v not in chosen]
            chosen |= set(extra[:size - len(chosen)])
        return chosen

    return {"uniform": uniform, "degree_stratified": degree_stratified,
            "disagreement_ego": disagreement_ego}[kind]


def greedy_on_subset(gbe, gfe, gme, nodes):
    """Greedy disagreement removal + exhaustive expansion restricted to
    `nodes` (a set of local indices). Returns the size found."""
    sub = lambda es: {(i, j) for i, j in es if i in nodes and j in nodes}
    be, fe, me = sub(gbe), sub(gfe), sub(gme)
    active = set(nodes)
    from collections import defaultdict
    for _ in range(10000):
        ab = {(i, j) for i, j in be if i in active and j in active}
        af = {(i, j) for i, j in fe if i in active and j in active}
        am = {(i, j) for i, j in me if i in active and j in active}
        dis = (ab ^ af) | (af ^ am) | (ab ^ am)
        if not dis:
            break
        sc = defaultdict(int)
        for u, v in dis:
            if u in active:
                sc[u] += 1
            if v in active:
                sc[v] += 1
        active.remove(max(sc, key=sc.get))
    # exhaustive expansion
    for nd in nodes:
        if nd in active:
            continue
        t = active | {nd}
        if ({(i, j) for i, j in be if i in t and j in t} ==
                {(i, j) for i, j in fe if i in t and j in t} ==
                {(i, j) for i, j in me if i in t and j in t}):
            active.add(nd)
    return len(active)


def exact_on_subset(gbe, gfe, gme, nodes):
    """Solve MCIS exactly via ILP (Max Independent Set on disagreement
    graph) with CBC. Returns (optimal_size, solve_seconds)."""
    nodes = sorted(nodes)
    sub = lambda es: {(i, j) for i, j in es if i in nodes and j in nodes}
    be, fe, me = sub(gbe), sub(gfe), sub(gme)
    # disagreement edges as undirected pairs
    union = be | fe | me
    dis_pairs = set()
    for (i, j) in union:
        if not ((i, j) in be and (i, j) in fe and (i, j) in me):
            dis_pairs.add((min(i, j), max(i, j)) if i != j else (i, j))
    prob = pulp.LpProblem("MCIS", pulp.LpMaximize)
    x = {v: pulp.LpVariable(f"x_{v}", cat="Binary") for v in nodes}
    prob += pulp.lpSum(x.values())
    for (i, j) in dis_pairs:
        if i == j:
            prob += x[i] == 0  # self-disagreement: node cannot be included
        else:
            prob += x[i] + x[j] <= 1
    t0 = time.time()
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    dt = time.time() - t0
    opt = int(round(pulp.value(prob.objective)))
    return opt, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--tiers", type=int, nargs="+", default=[20, 30, 40, 50])
    ap.add_argument("--per-tier", type=int, nargs="+", default=[15, 15, 10, 10])
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--sampler", default="uniform",
                    choices=["uniform", "degree_stratified", "disagreement_ego"],
                    help="subgraph sampling scheme (see make_sampler docstring)")
    ap.add_argument("--out", default=None,
                    help="output json filename (default ilp_validation.json)")
    args = ap.parse_args()

    data_dir = resolve_data_dir(args.data_dir)
    print(f"Data directory: {data_dir}")
    solver = build_solver(data_dir, n_seeds=1)
    solver._load()
    gbe, gfe, gme, ng = solver._gbe, solver._gfe, solver._gme, solver._ng
    rng = np.random.default_rng(args.seed)
    adj, _ = build_disagreement(set(gbe), set(gfe), set(gme), ng)
    sampler = make_sampler(adj, ng, args.sampler)
    print(f"Sampler: {args.sampler}")

    rows, all_gaps = [], []
    for size, count in zip(args.tiers, args.per_tier):
        g_ns, e_ns, times = [], [], []
        for _ in range(count):
            nodes = sampler(rng, size)
            gn = greedy_on_subset(gbe, gfe, gme, nodes)
            en, dt = exact_on_subset(gbe, gfe, gme, nodes)
            g_ns.append(gn)
            e_ns.append(en)
            times.append(dt)
            all_gaps.append((en - gn) / en if en > 0 else 0.0)
        g_ns, e_ns = np.array(g_ns), np.array(e_ns)
        gap = float((e_ns - g_ns).sum() / e_ns.sum()) if e_ns.sum() else 0.0
        rows.append({
            "size": size, "instances": count,
            "greedy_mean": float(g_ns.mean()), "greedy_sd": float(g_ns.std()),
            "exact_mean": float(e_ns.mean()), "exact_sd": float(e_ns.std()),
            "gap_pct": round(100 * gap, 2),
            "solve_time_s": [round(t, 4) for t in times],
        })
        print(f"  {size:>3d} nodes | {count:>2d} inst | "
              f"greedy {g_ns.mean():.1f}±{g_ns.std():.1f} | "
              f"exact {e_ns.mean():.1f}±{e_ns.std():.1f} | "
              f"gap {100*gap:.1f}% | t {min(times):.2f}-{max(times):.2f}s")

    all_gaps = np.array(all_gaps)
    summary = {
        "sampler": args.sampler,
        "tiers": rows,
        "overall": {
            "n_instances": int(sum(args.per_tier)),
            "mean_gap_pct": round(100 * float(all_gaps.mean()), 2),
            "max_gap_pct": round(100 * float(all_gaps.max()), 2),
            "frac_within_2pct": float((all_gaps <= 0.02).mean()),
        },
    }
    print(f"\n  Overall: mean gap {summary['overall']['mean_gap_pct']}%, "
          f"max gap {summary['overall']['max_gap_pct']}%, "
          f"{summary['overall']['frac_within_2pct']*100:.0f}% within 2%")

    (REPO / "results").mkdir(exist_ok=True)
    out = REPO / "results" / (args.out or "ilp_validation.json")
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Wrote {out}")


if __name__ == "__main__":
    main()
