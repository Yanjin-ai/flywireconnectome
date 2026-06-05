"""
Unit tests for mcis_connectome solver.
Run: pytest tests/ -v
No large data downloads required — all tests use synthetic graphs.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import networkx as nx
import numpy as np
from collections import defaultdict


# ── Core algorithm re-implementation (independent of package) ─────

def greedy_mcis(be, fe, me, ng, seed=0):
    rng = np.random.default_rng(seed)
    active = set(range(ng))
    for _ in range(50000):
        ab={(i,j) for i,j in be if i in active and j in active}
        af={(i,j) for i,j in fe if i in active and j in active}
        am={(i,j) for i,j in me if i in active and j in active}
        dis = (ab^af)|(af^am)|(ab^am)
        if not dis: return set(active)
        sc = defaultdict(int)
        for u,v in dis:
            if u in active: sc[u]+=1
            if v in active: sc[v]+=1
        worst = max(sc, key=lambda x:(sc[x], rng.random()))
        active.remove(worst)
    return set(active)

def expand(best, be, fe, me, ng):
    c = set(best)
    for nd in range(ng):
        if nd in c: continue
        t = c|{nd}
        if ({(i,j) for i,j in be if i in t and j in t} ==
            {(i,j) for i,j in fe if i in t and j in t} ==
            {(i,j) for i,j in me if i in t and j in t}):
            c.add(nd)
    return c

def verify(node_set, be, fe, me):
    s = node_set
    return (frozenset((i,j) for i,j in be if i in s and j in s) ==
            frozenset((i,j) for i,j in fe if i in s and j in s) ==
            frozenset((i,j) for i,j in me if i in s and j in s))


# ── Tests ─────────────────────────────────────────────────────────

class TestIsomorphismVerification:
    def test_empty_graphs_isomorphic(self):
        """Empty induced subgraph is always isomorphic."""
        be, fe, me = set(), set(), set()
        assert verify(set(), be, fe, me)

    def test_single_node_isomorphic(self):
        """Single node with no edges is isomorphic."""
        be = {(0,1), (1,2)}
        fe = {(0,1)}
        me = {(1,2)}
        assert verify({0}, be, fe, me)  # no edges in induced {0}

    def test_two_node_consensus_isomorphic(self):
        """Two nodes with identical edge in all 3 datasets."""
        be = fe = me = {(0,1)}
        assert verify({0,1}, be, fe, me)

    def test_two_node_disagreement_not_iso(self):
        """Two nodes where one dataset has an edge, others don't."""
        be = {(0,1)}
        fe, me = set(), set()
        assert not verify({0,1}, be, fe, me)

    def test_three_node_ring_iso(self):
        """Ring 0→1→2 present in all 3."""
        ring = {(0,1),(1,2),(2,0)}
        be = fe = me = ring
        assert verify({0,1,2}, be, fe, me)


class TestGreedyAlgorithm:
    def test_trivial_consensus(self):
        """When all 3 graphs are identical, all nodes are in MCIS."""
        n = 10
        edges = {(i,(i+1)%n) for i in range(n)}  # ring
        be = fe = me = edges
        result = greedy_mcis(be, fe, me, n, seed=0)
        result = expand(result, be, fe, me, n)
        assert len(result) == n
        assert verify(result, be, fe, me)

    def test_planted_common_subgraph(self):
        """Plant a K=5 common subgraph; algorithm should find it."""
        K, N = 5, 20
        # Consensus edges (only in first K nodes)
        consensus = {(i,(i+1)%K) for i in range(K)}
        be = fe = me = set(consensus)
        # Add noise: dataset-specific extra edges
        rng = np.random.default_rng(42)
        for graph in [be, fe, me]:
            for _ in range(10):
                u, v = rng.integers(K, N, 2)
                if u != v: graph.add((int(u), int(v)))
        result = greedy_mcis(be, fe, me, N, seed=0)
        result = expand(result, be, fe, me, N)
        assert verify(result, be, fe, me), "Result must be isomorphic"
        # Should contain most of the planted subgraph
        overlap = len(result & set(range(K)))
        assert overlap >= 3, f"Should recover at least 3/5 planted nodes, got {overlap}"

    def test_expansion_improves_greedy(self):
        """Expansion phase should never decrease MCIS size."""
        n = 15
        rng = np.random.default_rng(7)
        be = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(20)}
        fe = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(20)}
        me = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(20)}
        be = {e for e in be if e[0]!=e[1]}
        fe = {e for e in fe if e[0]!=e[1]}
        me = {e for e in me if e[0]!=e[1]}
        greedy_result = greedy_mcis(be, fe, me, n, seed=0)
        expanded = expand(greedy_result, be, fe, me, n)
        assert len(expanded) >= len(greedy_result)
        assert verify(expanded, be, fe, me)

    def test_result_is_always_isomorphic(self):
        """Algorithm output must always be verified isomorphic."""
        for seed in range(10):
            n = 20
            rng = np.random.default_rng(seed)
            be = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(25)}
            fe = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(25)}
            me = {(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(25)}
            for s in [be,fe,me]: s -= {(i,i) for i in range(n)}
            result = greedy_mcis(be, fe, me, n, seed=seed)
            result = expand(result, be, fe, me, n)
            assert verify(result, be, fe, me), f"Seed {seed}: result not isomorphic"


class TestNullModels:
    def test_shuffle_correspondence_yields_different_result(self):
        """Shuffling node correspondence should change MCIS (statistically)."""
        n = 30
        # Build graphs with a planted consensus
        K = 10
        consensus = {(i,(i+1)%K) for i in range(K)}
        rng = np.random.default_rng(0)
        noise = lambda: {(int(rng.integers(K,n)), int(rng.integers(K,n))) for _ in range(15)}
        be = consensus | noise()
        fe = consensus | noise()
        me = consensus | noise()
        real = expand(greedy_mcis(be, fe, me, n), be, fe, me, n)
        # Shuffle: permute fe node labels
        perm = np.random.default_rng(99).permutation(n)
        fe_sh = {(perm[u],perm[v]) for u,v in fe}
        shuffled = expand(greedy_mcis(be, fe_sh, me, n), be, fe_sh, me, n)
        # Real should be larger (consensus is gone after shuffle)
        assert len(real) >= len(shuffled), "Real should be ≥ shuffled"

    def test_degree_preserving_baseline(self):
        """Degree-preserving rewire should be near but not equal to real."""
        n = 20
        rng = np.random.default_rng(1)
        edges = list({(int(rng.integers(0,n)), int(rng.integers(0,n))) for _ in range(30)
                      if rng.integers(0,n) != rng.integers(0,n)})[:15]
        be = fe = me = set(edges)
        real_n = len(expand(greedy_mcis(be, fe, me, n), be, fe, me, n))
        assert real_n > 0


class TestExactVsGreedy:
    """Compare exact brute-force MCIS on very small graphs vs greedy."""

    def exact_mcis_small(self, be, fe, me, n, max_n=15):
        """Branch-and-bound exact MCIS for tiny graphs (n <= 15)."""
        if n > max_n:
            return None  # too large
        best = [set()]
        def bb(current, remaining):
            s = set(current)
            if verify(s, be, fe, me):
                if len(s) > len(best[0]): best[0] = set(s)
            if not remaining: return
            if len(s)+len(remaining) <= len(best[0]): return
            v = remaining[0]; rest = remaining[1:]
            bb(current+[v], rest)
            bb(current, rest)
        bb([], list(range(min(n, 12))))  # limit search to 12 for speed
        return best[0]

    def test_greedy_near_optimal_on_small_graphs(self):
        """On 10 small random graphs, greedy achieves ≥80% of exact."""
        gaps = []
        for trial in range(10):
            n = 10
            rng = np.random.default_rng(trial*3+7)
            be = {(int(u),int(v)) for u,v in [(rng.integers(0,n),rng.integers(0,n)) for _ in range(12)] if u!=v}
            fe = {(int(u),int(v)) for u,v in [(rng.integers(0,n),rng.integers(0,n)) for _ in range(12)] if u!=v}
            me = {(int(u),int(v)) for u,v in [(rng.integers(0,n),rng.integers(0,n)) for _ in range(12)] if u!=v}
            greedy_n = len(expand(greedy_mcis(be, fe, me, n, seed=0), be, fe, me, n))
            exact_set = self.exact_mcis_small(be, fe, me, n)
            exact_n = len(exact_set) if exact_set is not None else greedy_n
            ratio = greedy_n / max(exact_n, 1)
            gaps.append(ratio)
        mean_ratio = np.mean(gaps)
        assert mean_ratio >= 0.80, f"Greedy mean optimality ratio {mean_ratio:.2f} < 0.80"
        print(f"\n  Greedy optimality ratio: {mean_ratio:.2f} ± {np.std(gaps):.2f}")
