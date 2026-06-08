"""
Full-graph exact MIS certificate for N (upgrades "near-optimal" -> bounds).
==========================================================================

science.md §9 states that an exact ILP on the full 987-node instance is
"intractable", so N=105 is reported only as a greedy local optimum validated
on subgraphs. This script tries to do better: it solves (or bounds) the
Maximum Independent Set on the FULL disagreement graph D, turning the headline
from "near-optimal" into a certified interval  LB <= N <= UB.

Exactness levers used:
  1. Separability — MIS decomposes over connected components of D; every node
     with NO disagreement edge is in the MIS for free.
  2. Forced-out nodes — a node with a self-disagreement (an edge present in
     some but not all datasets as a self-loop) can never be in S; drop it.
  3. Warm start — feed the greedy 105-node solution to CBC as an incumbent
     lower bound so branch-and-bound only has to prove/raise it.
  4. Time-limited CBC — even without closing the gap, CBC returns a valid
     dual UPPER BOUND, so we always emit a rigorous certificate.

If CBC finds a feasible set LARGER than 105, that is the new true N and all
docs must be updated to it (user directive: align to true values).

Run:  MCIS_DATA_DIR=/path/to/data python src/exact_full_mis.py --max-seconds 600
Output: results/exact_full_mis.json
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import networkx as nx
import pulp

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir  # noqa
from run_analysis import build_solver  # noqa
from incremental_mcis import build_disagreement  # noqa


def greedy_warmstart(adj, nodes):
    """Max-degree greedy removal on the induced disagreement subgraph -> a
    feasible independent set (the production rule), used as CBC incumbent."""
    active = set(nodes)
    while True:
        d = {v: sum(1 for w in adj.get(v, ()) if w in active) for v in active}
        if not d or max(d.values()) == 0:
            break
        active.discard(max(d, key=d.get))
    return active


def solve_component(comp, adj, max_seconds, warm):
    """Exact/bounded MIS on one component. Returns (lb, ub, proven, chosen)."""
    comp = sorted(comp)
    if len(comp) == 1:
        return 1, 1, True, set(comp)
    prob = pulp.LpProblem("MIS", pulp.LpMaximize)
    x = {v: pulp.LpVariable(f"x_{v}", cat="Binary") for v in comp}
    prob += pulp.lpSum(x.values())
    seen = set()
    for v in comp:
        for w in adj.get(v, ()):
            if w in x and (w, v) not in seen:
                prob += x[v] + x[w] <= 1
                seen.add((v, w))
    # NOTE: we deliberately do NOT use CBC warmStart/setInitialValue here.
    # In this pulp+CBC build, warmStart treats the initial assignment as variable
    # *bounds*, fixing the 0-valued vars and returning a spurious "proven optimal"
    # BELOW a known feasible solution (it reported 82 while a feasible IS of 100
    # was supplied). Clean solve + time limit is correct; CBC still returns a
    # valid dual bound if the limit is hit.
    cmd = pulp.PULP_CBC_CMD(msg=0, timeLimit=max_seconds)
    prob.solve(cmd)
    chosen = {v for v in comp if x[v].value() and x[v].value() > 0.5}
    lb = len(chosen)
    # CBC exposes the best dual bound; fall back to LP relaxation otherwise
    status = pulp.LpStatus[prob.status]
    proven = status == "Optimal"
    ub = lb if proven else _lp_upper_bound(comp, adj)
    return lb, ub, proven, chosen


def _lp_upper_bound(comp, adj):
    """LP-relaxation upper bound on MIS for a component (valid dual bound)."""
    prob = pulp.LpProblem("MIS_LP", pulp.LpMaximize)
    x = {v: pulp.LpVariable(f"y_{v}", lowBound=0, upBound=1) for v in comp}
    prob += pulp.lpSum(x.values())
    seen = set()
    for v in comp:
        for w in adj.get(v, ()):
            if w in x and (w, v) not in seen:
                prob += x[v] + x[w] <= 1
                seen.add((v, w))
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    return int(np.floor(pulp.value(prob.objective) + 1e-6))


def lovasz_theta(D, **scs_kw):
    """Lovász theta of D — a tighter upper bound on the independence number
    than clique cover (alpha(D) <= theta(D)). Requires cvxpy + an SDP solver
    (SCS); returns None if unavailable. NOTE: this is a *numerical* SDP bound
    (SCS tolerance ~1e-3), so we report alpha <= floor(theta + 1e-3); the
    clique-cover bound remains the purely combinatorial, solver-free guarantee."""
    try:
        import cvxpy as cp
    except Exception:
        return None
    nodes = list(D.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    rows, cols = [], []
    for u, v in D.edges():
        if u != v:
            rows.append(idx[u])
            cols.append(idx[v])
    X = cp.Variable((n, n), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 1, X[rows, cols] == 0]
    prob = cp.Problem(cp.Maximize(cp.sum(X)), cons)
    prob.solve(solver=cp.SCS, **scs_kw)
    return float(prob.value)


def _one_clique_cover(D, order):
    """Single greedy clique cover given a vertex visitation order -> #cliques."""
    uncovered = set(D.nodes())
    pos = {v: i for i, v in enumerate(order)}
    cliques = 0
    for seed in order:
        if seed not in uncovered:
            continue
        clique = {seed}
        cand = sorted((u for u in D.neighbors(seed) if u in uncovered),
                      key=lambda v: pos[v])
        for u in cand:
            if all(D.has_edge(u, c) for c in clique):
                clique.add(u)
        uncovered -= clique
        cliques += 1
    return cliques


def clique_cover_ub(D, restarts=40, seed=0):
    """Valid UPPER bound on the independence number via greedy clique cover:
    any independent set hits each clique at most once, so MIS <= #cliques in
    ANY clique cover. We take the MINIMUM over several vertex orderings
    (degree-descending + randomised restarts) to tighten the bound."""
    rng = np.random.default_rng(seed)
    nodes = list(D.nodes())
    best = _one_clique_cover(D, sorted(nodes, key=lambda v: D.degree(v),
                                       reverse=True))
    for _ in range(restarts):
        order = nodes.copy()
        rng.shuffle(order)
        best = min(best, _one_clique_cover(D, order))
    return best


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--seeds", type=int, default=100,
                    help="greedy multi-start seeds for the certified lower bound")
    ap.add_argument("--cbc-seconds", type=int, default=0,
                    help="optional clean-CBC cross-check budget (0 = skip; the "
                         "full 987-node MIS does not close in reasonable time)")
    ap.add_argument("--theta", action="store_true",
                    help="also compute the Lovász theta SDP bound (needs cvxpy"
                         "+SCS; tighter than clique cover but numerical)")
    args = ap.parse_args()

    from run_analysis import seed_distribution  # noqa
    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    ng = solver._ng
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    adj, forced_out = build_disagreement(be, fe, me, ng)

    candidates = [v for v in range(ng) if v not in forced_out]
    n_edges = sum(len(adj.get(v, ())) for v in candidates) // 2
    dens = n_edges / (len(candidates) * (len(candidates) - 1) / 2)
    print(f"  Disagreement graph D: {len(candidates)} candidate nodes "
          f"({len(forced_out)} forced out), {n_edges} edges, density {dens:.3f}")

    D = nx.Graph()
    D.add_nodes_from(candidates)
    for v in candidates:
        for w in adj.get(v, ()):
            if w not in forced_out and v < w:
                D.add_edge(v, w)

    # ---- certified LOWER bound: greedy multi-start MCIS, verified independent
    sizes, best = seed_distribution(solver, args.seeds)
    lb = len(best)
    internal = sum(1 for v in best for w in adj.get(v, ()) if w in best) // 2
    assert internal == 0, f"greedy LB not independent! {internal} internal edges"
    print(f"  Certified LOWER bound (greedy multi-start, verified independent): "
          f"N >= {lb}  (0 internal disagreement edges)")

    # ---- UPPER bound: greedy clique cover (MIS <= #cliques) -----------------
    t0 = time.time()
    ub_clique = clique_cover_ub(D)
    print(f"  Clique-cover UPPER bound: N <= {ub_clique}  ({time.time()-t0:.1f}s)")
    ub = ub_clique

    # ---- optional tighter UPPER bound: Lovász theta (numerical SDP) ---------
    theta = None
    if args.theta:
        t0 = time.time()
        th = lovasz_theta(D)
        if th is not None:
            theta = {"value": round(th, 2), "ub_floor": int(np.floor(th + 1e-3)),
                     "seconds": round(time.time() - t0, 1)}
            print(f"  Lovász theta: {th:.2f} -> N <= {theta['ub_floor']} "
                  f"(numerical SDP; {theta['seconds']}s)")
            ub = min(ub, theta["ub_floor"])
        else:
            print("  Lovász theta: cvxpy/SCS not available, skipped")

    # ---- optional clean-CBC cross-check (no warmStart) ---------------------
    cbc = None
    if args.cbc_seconds > 0:
        comp = max(nx.connected_components(D), key=len)
        c_lb, c_ub, proven, _ = solve_component(sorted(comp), adj,
                                                args.cbc_seconds, set())
        cbc = {"lb": c_lb, "ub": c_ub, "proven_optimal": proven,
               "seconds": args.cbc_seconds}
        print(f"  Clean-CBC ({args.cbc_seconds}s): MIS in [{c_lb}, {c_ub}] "
              f"{'PROVEN' if proven else '(time-limited)'}")

    print(f"\n  CERTIFICATE:  {lb} <= N <= {ub}")
    prev = 105
    if lb > prev:
        print(f"  *** Greedy N={lb} EXCEEDS previously reported {prev} — "
              f"update docs to true value ***")
    else:
        print(f"  N={lb} is a verified-feasible MCIS; exact optimum is bounded "
              f"in [{lb}, {ub}] (full-graph MIS does not close exactly).")

    out = {
        "candidate_nodes": len(candidates),
        "forced_out": len(forced_out),
        "disagreement_edges": n_edges,
        "density": round(dens, 4),
        "lower_bound": {"N": lb, "method": "greedy multi-start + expansion",
                        "verified_independent": True,
                        "seed_max": int(sizes.max()),
                        "seed_mean": float(sizes.mean())},
        "upper_bound": {"N": ub,
                        "clique_cover": ub_clique,
                        "lovasz_theta": theta,
                        "method": ("Lovász theta (numerical) ∧ clique cover"
                                   if theta else "greedy clique cover")},
        "cbc_crosscheck": cbc,
        "certificate": f"{lb} <= N <= {ub}",
        "previous_reported_N": prev,
        "note": ("Full-graph exact MIS does not close in reasonable time (one "
                 "dense 987-node component, 45k constraints). The earlier "
                 "'82 proven' was a pulp/CBC warmStart bug — warmStart fixed "
                 "the 0-vars, yielding a spurious optimum below a known feasible "
                 "solution. N=105 is verified feasible. Clique cover is the "
                 "solver-free combinatorial UB; Lovász theta (if computed) is a "
                 "tighter numerical SDP UB (SCS tol ~1e-3)."),
    }
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "exact_full_mis.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"  Wrote results/exact_full_mis.json")


if __name__ == "__main__":
    main()
