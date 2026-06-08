"""
Maximum WEAKLY-CONNECTED Common Induced Subgraph (the connectivity requirement).
================================================================================

The plain MCIS (= Maximum Independent Set on the disagreement graph) is dominated
by neurons with ZERO conserved edges — they are "trivially isomorphic" (no edges
to disagree about). On the real data the N=109 MCIS is 97 disconnected pieces, 87
of them isolated nodes. That is why the degree-sequence null explained ~93% of N
(§4.2): the node count was inflated by edge-less neurons.

The connectivity requirement removes exactly that artifact. We now demand the
discovered structure be **weakly connected via the conserved (consensus) edges**.
The problem becomes the **Maximum Connected Common Induced Subgraph**:

    maximise |S|  s.t.  (1) no disagreement edge within S  (independent in D), and
                        (2) the consensus edges within S connect S (connected in C).

This is NP-hard (it is a connected variant of MIS), but the optimum is small, so
a connected-growth multi-start search finds it reliably. Every candidate is
verified: weakly connected AND zero internal disagreement.

Run:  MCIS_DATA_DIR=/path/to/data python src/connected_mcis.py --restarts 8000
Output: results/connected_mcis.json, results/connected_circuit.csv
"""
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir  # noqa
from run_analysis import build_solver  # noqa
from incremental_mcis import build_disagreement  # noqa


def consensus_adj(consensus):
    C = defaultdict(set)
    for (i, j) in consensus:
        if i != j:
            C[i].add(j)
            C[j].add(i)
    return C


def grow(seed, Cadj, Dadj, forced, rng):
    """Greedily grow a connected, disagreement-free set from `seed`. Candidate
    selection prefers nodes that keep the most future-compatible options."""
    S = {seed}
    while True:
        cand = set()
        for u in S:
            for v in Cadj[u]:
                if v in S or v in forced:
                    continue
                if any(v in Dadj[w] for w in S):
                    continue
                cand.add(v)
        if not cand:
            break
        # pick a candidate that adds the most compatible consensus-neighbours
        # (look-ahead), randomised among ties to diversify restarts
        scored = []
        for v in cand:
            opts = sum(1 for w in Cadj[v]
                       if w not in S and w not in forced
                       and not any(w in Dadj[x] for x in S))
            scored.append((opts + rng.random(), v))
        S.add(max(scored)[1])
    return S


def regrow(S0, Cadj, Dadj, forced, rng):
    """Grow a connected disagreement-free set starting from an existing
    connected core S0 (used by iterated local search)."""
    S = set(S0)
    while True:
        cand = set()
        for u in S:
            for v in Cadj[u]:
                if v in S or v in forced:
                    continue
                if any(v in Dadj[w] for w in S):
                    continue
                cand.add(v)
        if not cand:
            break
        scored = [(sum(1 for w in Cadj[v] if w not in S and w not in forced
                       and not any(w in Dadj[x] for x in S)) + rng.random(), v)
                  for v in cand]
        S.add(max(scored)[1])
    return S


def connected_core(S, consensus, keep, rng):
    """Return a connected sub-core of S of size ~keep (BFS from a random node)."""
    adj = defaultdict(set)
    Sset = set(S)
    for (i, j) in consensus:
        if i in Sset and j in Sset:
            adj[i].add(j); adj[j].add(i)
    start = rng.choice(list(S))
    seen, frontier = {int(start)}, [int(start)]
    while frontier and len(seen) < keep:
        u = frontier.pop(0)
        nbrs = list(adj[u]); rng.shuffle(nbrs)
        for v in nbrs:
            if v not in seen:
                seen.add(v); frontier.append(v)
                if len(seen) >= keep:
                    break
    return seen


def verify(S, consensus, Dadj):
    G = nx.Graph()
    G.add_nodes_from(S)
    for (i, j) in consensus:
        if i in S and j in S:
            G.add_edge(i, j)
    connected = nx.is_connected(G) if len(S) > 1 else True
    dis = sum(1 for u in S for v in Dadj[u] if v in S) // 2
    return connected, dis, G.number_of_edges()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--restarts", type=int, default=8000)
    ap.add_argument("--ils", type=int, default=20000,
                    help="iterated-local-search iterations after multi-start")
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    be, fe, me = set(solver._gbe), set(solver._gfe), set(solver._gme)
    ng = solver._ng
    consensus = be & fe & me
    Dadj, forced = build_disagreement(be, fe, me, ng)
    Cadj = consensus_adj(consensus)
    seeds = [i for i in range(ng) if i in Cadj and i not in forced]
    print(f"  Consensus graph: {len(consensus)} edges, {len(seeds)} candidate "
          f"seed nodes (have a conserved edge)")

    rng = np.random.default_rng(0)
    best = set()
    for r in range(args.restarts):
        S = grow(int(rng.choice(seeds)), Cadj, Dadj, forced, rng)
        if len(S) > len(best):
            best = set(S)
            print(f"    restart {r:5d}: NEW BEST connected N = {len(best)}")

    # ---- iterated local search: perturb (keep a connected core) + regrow -----
    for it in range(args.ils):
        keep = max(2, len(best) - int(rng.integers(1, 5)))
        core = connected_core(best, consensus, keep, rng)
        S = regrow(core, Cadj, Dadj, forced, rng)
        if len(S) > len(best):
            best = set(S)
            print(f"    ILS {it:5d}: NEW BEST connected N = {len(best)}")

    connected, dis, n_edges = verify(best, consensus, Dadj)
    assert connected and dis == 0, "invalid connected solution!"
    print(f"\n  Maximum connected common induced subgraph: N = {len(best)}, "
          f"{n_edges} conserved edges, weakly connected, 0 disagreement")

    # map to cell types
    gl, trip = solver._gl, solver._triples
    rows = []
    for loc in sorted(best):
        r = trip.iloc[gl[loc]]
        rows.append({"BANC": str(r["root_626"]), "cell_type": r.get("cell_type"),
                     "super_class": r.get("super_class"),
                     "neurotransmitter": r.get("neurotransmitter_predicted")})
    import pandas as pd
    os.makedirs(results_dir(), exist_ok=True)
    pd.DataFrame(rows).to_csv(results_dir() + "connected_circuit.csv", index=False)

    # write the deliverable network.csv + network_enriched.csv (the 27-node circuit)
    from run_analysis import write_enriched, REPO
    trip = solver._triples.iloc[[gl[i] for i in sorted(best)]].reset_index(drop=True)
    pd.DataFrame({"BANC": trip["root_626"].astype(str),
                  "FAFB": trip["fafb_match"].astype(str),
                  "MANC": trip["manc_match"].astype(str)}).to_csv(
        str(REPO / "network.csv"), index=False)
    write_enriched(trip, REPO / "network_enriched.csv")
    print(f"  Wrote network.csv + network_enriched.csv (the {len(best)}-node "
          f"connected deliverable)")

    out = {
        "N_connected": len(best),
        "conserved_edges": n_edges,
        "weakly_connected": bool(connected),
        "internal_disagreement": int(dis),
        "vs_unconstrained_MCIS": {"N": 109, "components": 97,
                                  "isolated_nodes": 87,
                                  "note": "plain MCIS is dominated by edge-less "
                                          "neurons; connectivity removes them"},
        "restarts": args.restarts,
    }
    with open(results_dir() + "connected_mcis.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"  Wrote results/connected_mcis.json, results/connected_circuit.csv")


if __name__ == "__main__":
    main()
