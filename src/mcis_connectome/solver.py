"""
MCISSolver: Maximum Common Induced Subgraph across connectomes.

Algorithm: Greedy Disagreement Removal + Exhaustive Expansion
Reference: See science.md §3.3 for full specification.

Complexity: O(N · D) per iteration where D = disagreement edges.
Guarantee: Local optimum; expansion is exhaustive.
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
        self.node_indices  = node_indices
        self.triples       = triples_df.iloc[sorted(node_indices)].copy()
        self.n             = len(node_indices)
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
            active.remove(max(sc, key=sc.get))
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

    def solve(self, n_seeds: Optional[int] = None) -> MCISResult:
        """
        Run MCIS solver with multi-seed greedy + exhaustive expansion.

        Returns
        -------
        MCISResult with the largest found common induced subgraph.
        """
        self._load()
        n_seeds = n_seeds or self.n_seeds
        best_set, best_n = set(), 0
        t0 = time.time()

        for seed in range(n_seeds):
            res, _ = self._greedy(self._gbe, self._gfe, self._gme, self._ng, seed)
            exp    = self._expand(res, self._gbe, self._gfe, self._gme, self._ng)
            if len(exp) > best_n:
                best_n, best_set = len(exp), exp
                self._log(f"  seed {seed:3d}: NEW BEST N={best_n}")

        self._log(f"  Solved in {time.time()-t0:.1f}s | best N={best_n}")

        global_idx = [self._gl[i] for i in sorted(best_set)]
        edge_sets  = [
            {(i,j) for i,j in self._gbe if i in best_set and j in best_set},
            {(i,j) for i,j in self._gfe if i in best_set and j in best_set},
            {(i,j) for i,j in self._gme if i in best_set and j in best_set},
        ]
        return MCISResult(best_set, self._triples.iloc[global_idx].reset_index(drop=True),
                          edge_sets, list(self.edge_lists.keys()))
