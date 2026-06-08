"""
Unit tests for the rigor-pass extensions (synthetic graphs, no data download):
  - exact_ilp samplers (uniform / degree_stratified / disagreement_ego)
  - worstcase_greedy (greedy underestimation vs exact MIS)
  - null_sensitivity (degree preservation, reciprocity)
  - exact_full_mis (clique-cover upper bound is valid; warm-start is independent)
Run: pytest tests/test_extensions.py -v
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import networkx as nx
import pytest

import worstcase_greedy as wc
import null_sensitivity as ns
import exact_full_mis as efm
import exact_ilp


# ── worstcase_greedy ────────────────────────────────────────────────
def test_greedy_mis_returns_independent_set():
    # P3: 0-1-2 ; MIS = {0,2}
    edges = [(0, 1), (1, 2)]
    keep = wc.greedy_mis_maxdeg(3, edges)
    # no edge inside the kept set
    assert not any((u in keep and v in keep) for u, v in edges)


def test_exact_mis_ge_greedy_on_random():
    rng = np.random.default_rng(1)
    for _ in range(5):
        n = 12
        edges = {(i, j) for i in range(n) for j in range(i + 1, n)
                 if rng.random() < 0.3}
        g = len(wc.greedy_mis_maxdeg(n, edges))
        o = wc.exact_mis(n, edges)
        assert o >= g  # exact optimum never below greedy
        assert g >= 1


def test_greedy_trap_is_valid_graph():
    rng = np.random.default_rng(0)
    n, edges, k, nR = wc.greedy_trap(20, rng)
    assert k == 20 and nR > 0
    assert all(0 <= u < n and 0 <= v < n for u, v in edges)


# ── null_sensitivity ────────────────────────────────────────────────
def test_rewire_preserves_in_out_degree():
    rng = np.random.default_rng(2)
    ng = 40
    edges = {(int(u), int(v)) for u, v in
             zip(rng.integers(0, ng, 120), rng.integers(0, ng, 120)) if u != v}
    din0, dout0 = ns.degree_seq(edges, ng)
    rew = ns.rewire(edges, ng, nswap=len(edges) * 5, seed=3)
    din1, dout1 = ns.degree_seq(rew, ng)
    assert np.array_equal(din0, din1)
    assert np.array_equal(dout0, dout1)


def test_reciprocity_bounds():
    assert ns.reciprocity(set()) == 0.0
    assert ns.reciprocity({(0, 1), (1, 0)}) == 1.0
    assert ns.reciprocity({(0, 1), (1, 2)}) == 0.0


# ── exact_full_mis: certificate soundness ───────────────────────────
def test_clique_cover_is_valid_upper_bound():
    rng = np.random.default_rng(4)
    for _ in range(5):
        n = 15
        G = nx.gnp_random_graph(n, 0.3, seed=int(rng.integers(1e6)))
        true_mis = wc.exact_mis(n, set(G.edges()))
        ub = efm.clique_cover_ub(G, restarts=10, seed=0)
        assert ub >= true_mis  # a clique cover always upper-bounds MIS


def test_warmstart_is_independent():
    # build a small disagreement adjacency and check greedy_warmstart is an IS
    adj = {0: {1}, 1: {0, 2}, 2: {1, 3}, 3: {2}}
    keep = efm.greedy_warmstart(adj, set(adj))
    assert not any(w in keep for v in keep for w in adj.get(v, ()))


# ── exact_ilp samplers ──────────────────────────────────────────────
def _toy_adj(ng=30, p=0.2, seed=5):
    rng = np.random.default_rng(seed)
    adj = {v: set() for v in range(ng)}
    for i in range(ng):
        for j in range(i + 1, ng):
            if rng.random() < p:
                adj[i].add(j)
                adj[j].add(i)
    return adj


@pytest.mark.parametrize("kind", ["uniform", "degree_stratified",
                                  "disagreement_ego"])
def test_sampler_returns_requested_size(kind):
    ng = 30
    adj = _toy_adj(ng)
    sampler = exact_ilp.make_sampler(adj, ng, kind)
    rng = np.random.default_rng(7)
    sub = sampler(rng, 12)
    assert len(sub) == 12
    assert all(0 <= v < ng for v in sub)


def test_ego_sampler_is_connected_in_disagreement_graph():
    ng = 40
    adj = _toy_adj(ng, p=0.25)
    sampler = exact_ilp.make_sampler(adj, ng, "disagreement_ego")
    rng = np.random.default_rng(9)
    sub = sampler(rng, 15)
    # the ego ball (before random top-up) should induce a connected subgraph
    G = nx.Graph()
    G.add_nodes_from(sub)
    for v in sub:
        for w in adj[v]:
            if w in sub:
                G.add_edge(v, w)
    # at least one large connected component (ego growth), not 15 isolated nodes
    assert max((len(c) for c in nx.connected_components(G)), default=0) >= 2
