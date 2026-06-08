"""
MCISSolver: Maximum Common Induced Subgraph across connectomes.

MCIS (with the node correspondence given) = Maximum Independent Set on the
disagreement graph D. Two solvers live here:

  * Production: GMIN (minimum-degree greedy SELECTION) + (1,2)-swap local search,
    multi-start. Finds N≈110 — strictly better than the baseline below, which
    systematically underestimates MIS (science.md §3.4, §3.6).
  * Baseline: max-degree disagreement removal + exhaustive expansion
    (`_greedy`/`_expand`) — kept as the documented vertex-cover heuristic that
    GMIN beats (it returns only ~105 on this instance).

Reference: science.md §3.2–3.6.
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
from typing import Dict, Optional, List
import time


class MCISResult:
    """Results from an MCIS solve."""
    def __init__(self, node_indices, triples_df, edge_sets, dataset_names):
        # triples_df is expected to already contain exactly the selected rows,
        # aligned to sorted(node_indices). We do NOT re-index it here.
        self.node_indices  = sorted(node_indices)
        self.triples       = triples_df.reset_index(drop=True).copy()
        self.n             = len(self.triples)
        self.dataset_names = dataset_names
        # Verify isomorphism and count edges
        induced = [{(i,j) for i,j in es if i in node_indices and j in node_indices}
                   for es in edge_sets]
        self.is_isomorphic = all(e == induced[0] for e in induced)
        self.n_edges       = len(induced[0])

    def to_csv(self, path: str):
        cols = [c for c in self.triples.columns
                if c in ['root_626','fafb_match','manc_match']]
        out = self.triples[cols].copy()
        out.columns = self.dataset_names
        out.to_csv(path, index=False)
        print(f"Saved {self.n} rows to {path}")
        return out

    def summary(self) -> str:
        lines = [
            f"MCIS Result",
            f"  N neurons:        {self.n}",
            f"  N conserved edges:{self.n_edges}",
            f"  Isomorphic:       {self.is_isomorphic}",
            f"  Datasets:         {' × '.join(self.dataset_names)}",
        ]
        if 'super_class' in self.triples.columns:
            sc = self.triples['super_class'].value_counts()
            lines.append(f"  Composition:      {dict(sc)}")
        return '\n'.join(lines)

    def __repr__(self):
        return f"MCISResult(N={self.n}, edges={self.n_edges}, isomorphic={self.is_isomorphic})"


class MCISSolver:
    """
    Find the Maximum Common Induced Subgraph across 3+ connectome datasets.

    Parameters
    ----------
    edge_lists : dict[str, str]
        Mapping of dataset_name → path to edge list CSV.
        CSV must have columns: 'source neuron id', 'target neuron id'
        (or any two-column CSV; first = source, second = target).

    triplets_path : str
        Path to BANC metadata feather file containing cross-dataset
        neuron correspondence columns (fafb_match, manc_match).
        Downloadable from:
        gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/

    n_seeds : int
        Number of random tie-breaking seeds for greedy algorithm.
        Higher = more likely to find global optimum. Default: 10.

    verbose : bool
        Print progress. Default: True.
    """

    def __init__(self,
                 edge_lists: Dict[str, str],
                 triplets_path: str,
                 n_seeds: int = 10,
                 verbose: bool = True):
        self.edge_lists    = edge_lists
        self.triplets_path = triplets_path
        self.n_seeds       = n_seeds
        self.verbose       = verbose
        self._loaded       = False

    def _log(self, msg):
        if self.verbose: print(msg)

    def _load(self):
        if self._loaded: return

        self._log("Loading BANC metadata...")
        bm = pd.read_feather(self.triplets_path)
        m  = bm[bm['fafb_match'].notna() & bm['manc_match'].notna()].copy()
        m  = m.drop_duplicates('root_626')
        for c in ['root_626','fafb_match','manc_match']:
            m[c] = m[c].astype(str)
        self._meta = m
        self._log(f"  {len(m):,} triplets with FAFB + MANC matches")

        def lf(fp, ns):
            d = pd.read_csv(fp, dtype=str)
            d.columns = ['source','target']
            d = d[d['source'].isin(ns) & d['target'].isin(ns)]
            return nx.from_pandas_edgelist(d,'source','target',
                                           create_using=nx.DiGraph())

        graphs = {}
        for name, path in self.edge_lists.items():
            self._log(f"Loading {name}...")
            if name == 'BANC':
                G = lf(path, set(m['root_626']))
            elif name == 'FAFB':
                G = lf(path, set(m['fafb_match']))
            elif name == 'MANC':
                G = lf(path, set(m['manc_match']))
            else:
                raise ValueError(f"Unknown dataset: {name}. Use BANC, FAFB, or MANC.")
            self._log(f"  {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
            graphs[name] = G

        ia3 = (set(graphs['BANC'].nodes()) &
               set(m[m['fafb_match'].isin(graphs['FAFB'].nodes())]['root_626']) &
               set(m[m['manc_match'].isin(graphs['MANC'].nodes())]['root_626']))
        self._triples = m[m['root_626'].isin(ia3)].reset_index(drop=True)
        self._log(f"  {len(self._triples):,} triplets in all 3 edge lists")

        bi = {r.root_626:   i for i,r in self._triples.iterrows()}
        fi = {r.fafb_match: i for i,r in self._triples.iterrows()}
        mi = {r.manc_match: i for i,r in self._triples.iterrows()}
        self._be = {(bi[u],bi[v]) for u,v in graphs['BANC'].edges() if u in bi and v in bi}
        self._fe = {(fi[u],fi[v]) for u,v in graphs['FAFB'].edges() if u in fi and v in fi}
        self._me = {(mi[u],mi[v]) for u,v in graphs['MANC'].edges() if u in mi and v in mi}

        ag = self._be & self._fe & self._me
        CG = nx.DiGraph(); CG.add_edges_from(ag)
        giant = max(nx.weakly_connected_components(CG), key=len)
        self._gl  = sorted(giant)
        self._ng  = len(self._gl)
        loc = {g:l for l,g in enumerate(self._gl)}
        self._gbe = {(loc[u],loc[v]) for u,v in self._be if u in loc and v in loc}
        self._gfe = {(loc[u],loc[v]) for u,v in self._fe if u in loc and v in loc}
        self._gme = {(loc[u],loc[v]) for u,v in self._me if u in loc and v in loc}
        self._log(f"  Giant consensus component: {self._ng} nodes")
        self._loaded = True

    @staticmethod
    def _greedy(gbe, gfe, gme, ng, seed=0):
        rng = np.random.default_rng(seed)
        active = set(range(ng))
        for _ in range(50000):
            ab = {e for e in gbe if e[0] in active and e[1] in active}
            af = {e for e in gfe if e[0] in active and e[1] in active}
            am = {e for e in gme if e[0] in active and e[1] in active}
            dis = (ab ^ af) | (af ^ am) | (ab ^ am)
            if not dis: return set(active), len(ab)
            sc = defaultdict(int)
            for u,v in dis:
                if u in active: sc[u] += 1
                if v in active: sc[v] += 1
            # randomised tie-breaking among the worst nodes (seed-dependent,
            # so multi-seed restarts explore genuinely different orderings)
            top = max(sc.values())
            worst = [n for n, s in sc.items() if s == top]
            active.remove(worst[int(rng.integers(len(worst)))])
        return set(active), -1

    @staticmethod
    def _expand(best, gbe, gfe, gme, ng):
        current = set(best)
        for nd in range(ng):
            if nd in current: continue
            t = current | {nd}
            if ({e for e in gbe if e[0] in t and e[1] in t} ==
                {e for e in gfe if e[0] in t and e[1] in t} ==
                {e for e in gme if e[0] in t and e[1] in t}):
                current.add(nd)
        return current

    # ── Stronger production solver: GMIN + (1,2)-swap on the disagreement graph
    # MCIS = Maximum Independent Set on the disagreement graph D (a pair of
    # neurons cannot coexist if their connection disagrees across connectomes).
    # The max-degree-removal greedy above is the vertex-cover heuristic (only
    # Θ(log n) for MIS) and *systematically underestimates* N on this instance
    # (science.md §3.4). GMIN (minimum-degree SELECTION; (Δ+2)/3 guarantee) plus
    # a (1,2)-swap local search finds a larger common subgraph — empirically
    # N≈110 vs the old 105. We keep _greedy/_expand as the documented baseline.

    @staticmethod
    def _disagreement_adj(gbe, gfe, gme, ng):
        """Undirected disagreement adjacency + nodes forced out (self-loop)."""
        union = gbe | gfe | gme
        consensus = gbe & gfe & gme
        adj = {v: set() for v in range(ng)}
        forced = set()
        for (i, j) in union:
            if (i, j) in consensus:
                continue
            if i == j:
                forced.add(i)
            else:
                adj[i].add(j)
                adj[j].add(i)
        return adj, forced

    @staticmethod
    def _gmin(adj, nodes, rng):
        """Minimum-degree greedy independent set on the disagreement graph."""
        active = set(nodes)
        deg = {v: sum(1 for w in adj[v] if w in active) for v in active}
        S = set()
        while active:
            mind = min(deg[v] for v in active)
            cands = [v for v in active if deg[v] == mind]
            v = cands[int(rng.integers(len(cands)))]
            S.add(v)
            gone = {v} | {w for w in adj[v] if w in active}
            active -= gone
            for g in gone:
                for w in adj[g]:
                    if w in active:
                        deg[w] -= 1
        return S

    @staticmethod
    def _two_swap(S, adj, nodes, rng, max_passes=6):
        """(1,2)-swap local search: drop one IS node, add two non-adjacent free
        nodes (net +1) where possible — escapes maximal-but-not-maximum sets."""
        S = set(S)
        node_list = list(nodes)
        for _ in range(max_passes):
            improved = False
            for v in list(S):
                Smv = S - {v}
                blocked = set()
                for u in Smv:
                    blocked |= adj[u]
                free = [w for w in node_list if w not in Smv and w not in blocked]
                if len(free) < 2:
                    continue
                rng.shuffle(free)
                cap = min(len(free), 120)
                found = False
                for i in range(cap):
                    a = free[i]; na = adj[a]
                    for j in range(i + 1, cap):
                        b = free[j]
                        if b not in na:
                            S = Smv | {a, b}; improved = found = True; break
                    if found:
                        break
                if found:
                    break
            if not improved:
                break
        return S

    @classmethod
    def _best_mcis(cls, gbe, gfe, gme, ng, restarts=2000, seed0=0):
        """Multi-start GMIN + (1,2)-swap → the largest common induced subgraph
        found. Returns (best_set, sizes_list)."""
        adj, forced = cls._disagreement_adj(gbe, gfe, gme, ng)
        nodes = [v for v in range(ng) if v not in forced]
        best, best_set, sizes = 0, set(), []
        for r in range(restarts):
            rng = np.random.default_rng(seed0 + r)
            S = cls._two_swap(cls._gmin(adj, nodes, rng), adj, nodes, rng)
            sizes.append(len(S))
            if len(S) > best:
                best, best_set = len(S), set(S)
        return best_set, sizes

    def solve(self, n_seeds: Optional[int] = None) -> MCISResult:
        """
        Run MCIS solver with multi-seed greedy + exhaustive expansion.

        Returns
        -------
        MCISResult with the largest found common induced subgraph.
        """
        self._load()
        n_seeds = n_seeds or self.n_seeds
        t0 = time.time()

        # Production solver: GMIN + (1,2)-swap multi-start (stronger than the
        # baseline max-degree greedy, which underestimates N; science.md §3.4).
        restarts = max(n_seeds * 30, 2000)
        best_set, _sizes = self._best_mcis(self._gbe, self._gfe, self._gme,
                                           self._ng, restarts=restarts)
        best_n = len(best_set)
        self._log(f"  Solved in {time.time()-t0:.1f}s | best N={best_n} "
                  f"(GMIN+2-swap, {restarts} restarts)")

        global_idx = [self._gl[i] for i in sorted(best_set)]
        edge_sets  = [
            {(i,j) for i,j in self._gbe if i in best_set and j in best_set},
            {(i,j) for i,j in self._gfe if i in best_set and j in best_set},
            {(i,j) for i,j in self._gme if i in best_set and j in best_set},
        ]
        return MCISResult(best_set, self._triples.iloc[global_idx].reset_index(drop=True),
                          edge_sets, list(self.edge_lists.keys()))
