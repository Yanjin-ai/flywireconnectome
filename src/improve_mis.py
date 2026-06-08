"""
Attack the 105–136 certificate gap from the LOWER-bound side (closes the
"NP-hardness used as a shield" gap).
=========================================================================

The production solver removes the MAX-degree disagreement node each step — the
classic vertex-cover heuristic, which is only a Θ(log n) approximation for MIS
(§3.4). Here we throw a *stronger* MIS attack at the full 987-node disagreement
graph to test whether N = 105 can be beaten:

  1. GMIN — minimum-degree greedy SELECTION (repeatedly take the min-degree
     vertex into the independent set, delete it + its neighbours). GMIN has a
     better guarantee for MIS ((Δ+2)/3) than max-degree removal, so it can find
     larger independent sets the production greedy misses.
  2. (1,2)-swap local search — remove one IS vertex and try to add two: the
     standard local move that escapes maximal-but-not-maximum independent sets.
  3. Many random restarts; report the best independent set found.

Plus a kernelization probe (degree-0/1 reductions) to see whether the instance
shrinks. If the best IS EXCEEDS 105, that is the new true N (and all docs must
follow). If it does NOT — after a strong, theory-backed search — that is solid
evidence the optimum is near 105, and the 105–136 gap is the loose UB, not a
missed larger circuit.

Run:  MCIS_DATA_DIR=/path/to/data python src/improve_mis.py --restarts 300
Output: results/improve_mis.json
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir  # noqa
from run_analysis import build_solver, seed_distribution  # noqa
from incremental_mcis import build_disagreement  # noqa


def gmin(adj, nodes, rng):
    """Minimum-degree greedy independent set."""
    active = set(nodes)
    deg = {v: sum(1 for w in adj[v] if w in active) for v in active}
    IS = set()
    while active:
        mind = min(deg[v] for v in active)
        cands = [v for v in active if deg[v] == mind]
        v = cands[int(rng.integers(len(cands)))]
        IS.add(v)
        # remove v and its neighbours from active, update degrees lazily
        gone = {v} | {w for w in adj[v] if w in active}
        active -= gone
        for g in gone:
            for w in adj[g]:
                if w in active:
                    deg[w] -= 1
    return IS


def is_independent(S, adj):
    return not any(w in S for v in S for w in adj[v])


def two_swap(IS, adj, nodes, rng, max_passes=6):
    """(1,2)-swap local search: drop one IS vertex, add two non-adjacent free
    vertices when possible (net +1). Returns an improved independent set."""
    IS = set(IS)
    node_list = list(nodes)
    for _ in range(max_passes):
        improved = False
        for v in list(IS):
            ISmv = IS - {v}
            blocked = set()
            for u in ISmv:
                blocked |= adj[u]
            free = [w for w in node_list if w not in ISmv and w not in blocked]
            if len(free) < 2:
                continue
            rng.shuffle(free)
            found = False
            # cap the pair search to keep it fast on dense free sets
            cap = min(len(free), 120)
            for i in range(cap):
                a = free[i]
                na = adj[a]
                for j in range(i + 1, cap):
                    b = free[j]
                    if b not in na:
                        IS = ISmv | {a, b}
                        improved = found = True
                        break
                if found:
                    break
            if found:
                break
        if not improved:
            break
    return IS


def kernelize_stats(adj, nodes):
    """Count how many vertices are removable by degree-0/1 reductions (a proxy
    for how much the dense instance kernelizes — expected small here)."""
    active = set(nodes)
    removed = 0
    changed = True
    while changed:
        changed = False
        for v in list(active):
            d = sum(1 for w in adj[v] if w in active)
            if d <= 1:  # degree-0 (free) or degree-1 (pendant: take v, drop nbr)
                active.discard(v)
                if d == 1:
                    nbr = next(w for w in adj[v] if w in active)
                    active.discard(nbr)
                removed += 1
                changed = True
    return removed, len(active)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--restarts", type=int, default=300)
    ap.add_argument("--seeds", type=int, default=40,
                    help="production-greedy multi-start for the baseline LB")
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    adj, forced = build_disagreement(set(solver._gbe), set(solver._gfe),
                                     set(solver._gme), ng)
    nodes = [v for v in range(ng) if v not in forced]
    for v in nodes:
        adj.setdefault(v, set())

    # baseline: the OLD max-degree-removal greedy + exhaustive expansion
    # (this is what underestimated N; GMIN below is now the production solver)
    base = 0
    for seed in range(args.seeds):
        res, _ = solver._greedy(solver._gbe, solver._gfe, solver._gme,
                                solver._ng, seed)
        exp = solver._expand(res, solver._gbe, solver._gfe, solver._gme,
                             solver._ng)
        base = max(base, len(exp))
    print(f"  Baseline max-degree-removal greedy ({args.seeds} seeds): "
          f"best N = {base}")

    # GMIN + 2-swap multi-start
    rng = np.random.default_rng(0)
    best, best_set = 0, set()
    t0 = time.time()
    gmin_best = 0
    for r in range(args.restarts):
        S = gmin(adj, nodes, rng)
        assert is_independent(S, adj)
        gmin_best = max(gmin_best, len(S))
        S = two_swap(S, adj, nodes, rng)
        assert is_independent(S, adj), "2-swap broke independence"
        if len(S) > best:
            best, best_set = len(S), set(S)
            print(f"    restart {r:3d}: NEW BEST IS = {best}  "
                  f"[{time.time()-t0:.0f}s]")
    print(f"  GMIN alone best: {gmin_best} | GMIN+2-swap best: {best} "
          f"({args.restarts} restarts, {time.time()-t0:.0f}s)")

    removed, kernel = kernelize_stats(adj, nodes)
    print(f"  Kernelization (deg-0/1): removed {removed}, residual core "
          f"{kernel} nodes (dense graph -> little reduction expected)")

    prev = 105
    overall = max(best, base)
    verdict = (f"best IS found = {overall}; "
               + ("EXCEEDS 105 — new true N, update docs"
                  if overall > prev else
                  f"does NOT beat {prev} despite a stronger (GMIN+2-swap) "
                  f"search — strong evidence the optimum is near {prev}; the "
                  f"105–136 gap is the loose Lovász-θ upper bound, not a missed "
                  f"larger circuit."))
    print(f"  VERDICT: {verdict}")

    out = {
        "production_greedy_best": int(base),
        "gmin_best": int(gmin_best),
        "gmin_plus_2swap_best": int(best),
        "overall_best_LB": int(overall),
        "restarts": args.restarts,
        "kernelization": {"removed_deg01": int(removed),
                          "residual_core_nodes": int(kernel)},
        "upper_bound_lovasz_theta": 136,
        "gap_after": f"{overall} <= N <= 136",
        "previous_reported_N": prev,
        "verdict": verdict,
    }
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "improve_mis.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  Wrote results/improve_mis.json")


if __name__ == "__main__":
    main()
