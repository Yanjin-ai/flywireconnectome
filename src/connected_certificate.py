"""
Optimality certificate for the connected MCIS (closes review critique P2).
==========================================================================

The connected deliverable (N=27) is a strong lower bound but had no upper bound.
Here we bracket it rigorously.

Key structural fact: a weakly-connected disagreement-free set is connected via
"clean" consensus edges (consensus edges whose endpoints do not disagree), so it
lives entirely inside ONE connected component of the *clean consensus graph*.
Therefore:

  N_connected  ≤  max over clean-components K of  α(D restricted to K)
               ≤  max over K of  ϑ(D|K)            (Lovász theta, an SDP UB)

The clean components are tiny except one giant component; we compute ϑ on the
giant one (cvxpy/SCS) for a rigorous upper bound far tighter than the trivial
α(D)=109. We ALSO solve a single-commodity-flow connectivity ILP on the giant
clean component (CBC, time-limited): its best integer solution is a *certified
feasible* connected circuit (a lower bound, possibly > 27 → would update N), and
if CBC proves optimality the connected MCIS is solved exactly.

Run:  MCIS_DATA_DIR=/path python src/connected_certificate.py --ilp-seconds 600
Output: results/connected_certificate.json
"""
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir  # noqa
from run_analysis import build_solver  # noqa
from incremental_mcis import build_disagreement  # noqa


def clean_consensus_graph(consensus, Dadj):
    """Consensus edges between disagreement-free pairs (edges a connected
    disagreement-free circuit is allowed to use)."""
    G = nx.Graph()
    for (i, j) in consensus:
        if i != j and j not in Dadj[i]:
            G.add_edge(i, j)
    return G


def theta_ub(nodes, Dadj):
    """Lovász theta of D restricted to `nodes` -> rigorous UB on alpha there."""
    try:
        import cvxpy as cp
    except Exception:
        return None
    nodes = sorted(nodes)
    idx = {v: k for k, v in enumerate(nodes)}
    n = len(nodes)
    rows, cols = [], []
    for v in nodes:
        for u in Dadj[v]:
            if u in idx:
                rows.append(idx[v]); cols.append(idx[u])
    X = cp.Variable((n, n), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 1]
    if rows:
        cons.append(X[rows, cols] == 0)
    prob = cp.Problem(cp.Maximize(cp.sum(X)), cons)
    prob.solve(solver=cp.SCS)
    return float(prob.value)


def scf_connectivity_ilp(comp, consensus, Dadj, seconds):
    """Single-commodity-flow ILP for the maximum weakly-connected disagreement-
    free induced subgraph within `comp`. Returns (lb_feasible, proven_optimal).
    A virtual source feeds exactly one root; each chosen node consumes 1 unit;
    flow only travels between chosen nodes -> selection must be connected."""
    import pulp
    comp = sorted(comp)
    cset = set(comp)
    # undirected clean edges within comp
    edges = set()
    for (i, j) in consensus:
        if i in cset and j in cset and i != j and j not in Dadj[i]:
            edges.add((min(i, j), max(i, j)))
    n = len(comp)
    prob = pulp.LpProblem("connected_mcis", pulp.LpMaximize)
    x = {v: pulp.LpVariable(f"x_{v}", cat="Binary") for v in comp}
    g = {v: pulp.LpVariable(f"g_{v}", cat="Binary") for v in comp}  # root
    # bidirected flow arcs
    arcs = []
    for (i, j) in edges:
        arcs += [(i, j), (j, i)]
    f = {a: pulp.LpVariable(f"f_{a[0]}_{a[1]}", lowBound=0) for a in arcs}
    fs = {v: pulp.LpVariable(f"fs_{v}", lowBound=0) for v in comp}  # source->root

    prob += pulp.lpSum(x.values())
    # disagreement independence
    seen = set()
    for v in comp:
        for u in Dadj[v]:
            if u in x and (u, v) not in seen:
                prob += x[v] + x[u] <= 1
                seen.add((v, u))
    # exactly one root, root must be chosen
    prob += pulp.lpSum(g.values()) == 1
    for v in comp:
        prob += g[v] <= x[v]
        prob += fs[v] <= n * g[v]
    # source sends total = number chosen, only via the root arc
    prob += pulp.lpSum(fs.values()) == pulp.lpSum(x.values())
    # flow conservation: each chosen node consumes exactly 1
    out_arcs = defaultdict(list)
    in_arcs = defaultdict(list)
    for a in arcs:
        out_arcs[a[0]].append(a)
        in_arcs[a[1]].append(a)
    for v in comp:
        prob += (fs[v] + pulp.lpSum(f[a] for a in in_arcs[v])
                 - pulp.lpSum(f[a] for a in out_arcs[v]) == x[v])
    # capacity: flow only between chosen nodes
    for a in arcs:
        prob += f[a] <= n * x[a[0]]
        prob += f[a] <= n * x[a[1]]

    prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=seconds))
    status = pulp.LpStatus[prob.status]
    lb = int(round(sum(1 for v in comp if x[v].value() and x[v].value() > 0.5)))
    return lb, status == "Optimal", n, len(edges)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--ilp-seconds", type=int, default=600)
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    consensus = be & fe & me
    Dadj, forced = build_disagreement(be, fe, me, solver._ng)

    clean = clean_consensus_graph(consensus, Dadj)
    comps = sorted(nx.connected_components(clean), key=len, reverse=True)
    giant = comps[0]
    print(f"  Clean-consensus components: sizes {[len(c) for c in comps[:6]]}; "
          f"giant = {len(giant)} nodes")

    # ---- rigorous UPPER bound: theta on the giant clean component -----------
    t0 = time.time()
    th = theta_ub(giant, Dadj)
    ub_theta = int(np.floor(th + 1e-3)) if th is not None else None
    print(f"  ϑ(D | giant clean component) = {th:.2f} -> N_connected ≤ {ub_theta}"
          f"  ({time.time()-t0:.0f}s)") if th else print("  (cvxpy unavailable)")

    # ---- connectivity ILP: certified feasible LB (+ exact if it closes) -----
    t0 = time.time()
    lb, proven, n_ilp, m_ilp = scf_connectivity_ilp(giant, consensus, Dadj,
                                                    args.ilp_seconds)
    print(f"  Connectivity ILP on {n_ilp}-node / {m_ilp}-edge clean component: "
          f"feasible N = {lb} {'(PROVEN OPTIMAL)' if proven else '(time-limited)'}"
          f"  ({time.time()-t0:.0f}s)")

    lb_final = max(27, lb)
    ub_final = ub_theta if ub_theta is not None else 109
    if proven:
        ub_final = lb  # exact
    print(f"\n  CONNECTED CERTIFICATE:  {lb_final} ≤ N_connected ≤ {ub_final}"
          f"  {'(EXACT)' if proven else ''}")
    if lb > 27:
        print(f"  *** ILP found a LARGER connected circuit ({lb} > 27) — update ***")

    out = {
        "heuristic_lb": 27,
        "clean_giant_component": len(giant),
        "theta_ub_giant": th, "theta_ub_floor": ub_theta,
        "ilp": {"feasible_lb": lb, "proven_optimal": proven,
                "nodes": n_ilp, "edges": m_ilp, "seconds": args.ilp_seconds},
        "certificate": {"lb": lb_final, "ub": ub_final,
                        "exact": bool(proven)},
        "note": ("A connected disagreement-free set lives in one clean-consensus "
                 "component, so ϑ on the giant clean component is a rigorous UB "
                 "(far tighter than the trivial α(D)=109). The single-commodity-"
                 "flow ILP gives a certified-feasible connected circuit and, if it "
                 "closes, the exact connected MCIS."),
    }
    os.makedirs(results_dir(), exist_ok=True)
    with open(results_dir() + "connected_certificate.json", "w") as fju:
        json.dump(out, fju, indent=2)
    print("  Wrote results/connected_certificate.json")


if __name__ == "__main__":
    main()
