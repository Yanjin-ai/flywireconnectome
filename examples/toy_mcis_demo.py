"""
Toy demonstration of the MCIS algorithm on small synthetic graphs.
No large datasets needed — runs in seconds.

This shows the core algorithm logic on a 20-node example where the
ground-truth MCIS is known, then scales to show how it behaves on
random graphs of increasing size.

Run: python examples/toy_mcis_demo.py
"""

import networkx as nx
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict


# ── Core MCIS algorithm (same as production code) ──────────────────

def greedy_mcis(G1, G2, G3, node_set, seed=0):
    """
    Find Maximum Common Induced Subgraph via greedy disagreement removal.

    Parameters
    ----------
    G1, G2, G3 : nx.DiGraph
        Three directed graphs over the SAME node set (bijection assumed).
    node_set : list
        Nodes to consider (indices into G1/G2/G3).
    seed : int
        Random seed for tie-breaking.

    Returns
    -------
    set of nodes forming the largest found common induced subgraph.
    """
    rng    = np.random.default_rng(seed)
    active = set(node_set)

    for _ in range(10000):
        # Disagreement = edge present in some but not all datasets
        e1 = {(u,v) for u,v in G1.edges() if u in active and v in active}
        e2 = {(u,v) for u,v in G2.edges() if u in active and v in active}
        e3 = {(u,v) for u,v in G3.edges() if u in active and v in active}
        dis = (e1^e2) | (e2^e3) | (e1^e3)

        if not dis:
            break  # Isomorphic — done

        # Score each node by number of disagreement edges it participates in
        score = defaultdict(int)
        for u, v in dis:
            if u in active: score[u] += 1
            if v in active: score[v] += 1

        # Remove node with most disagreements (randomised tie-breaking)
        worst = max(score, key=lambda x: (score[x], rng.random()))
        active.remove(worst)

    # Expansion: try adding back removed nodes
    removed = set(node_set) - active
    for nd in removed:
        test = active | {nd}
        e1t = {(u,v) for u,v in G1.edges() if u in test and v in test}
        e2t = {(u,v) for u,v in G2.edges() if u in test and v in test}
        e3t = {(u,v) for u,v in G3.edges() if u in test and v in test}
        if e1t == e2t == e3t:
            active.add(nd)

    return active


def verify_isomorphic(G1, G2, G3, node_set):
    """Check that induced subgraphs are identical across all three datasets."""
    s = node_set
    e1 = frozenset((u,v) for u,v in G1.edges() if u in s and v in s)
    e2 = frozenset((u,v) for u,v in G2.edges() if u in s and v in s)
    e3 = frozenset((u,v) for u,v in G3.edges() if u in s and v in s)
    return e1 == e2 == e3


# ── Example 1: Known ground truth ──────────────────────────────────

print("=" * 55)
print("Example 1: Synthetic graph with known ground truth")
print("=" * 55)

np.random.seed(42)
N = 20   # total nodes
K = 8    # ground-truth MCIS size

# Build "consensus" edges (present in all 3 datasets)
consensus_nodes = list(range(K))
consensus_edges = [(i, (i+1)%K) for i in range(K)]   # ring topology

# Build three graphs
G1, G2, G3 = nx.DiGraph(), nx.DiGraph(), nx.DiGraph()
for G in [G1, G2, G3]:
    G.add_nodes_from(range(N))
    G.add_edges_from(consensus_edges)  # shared edges

# Add noise: dataset-specific extra edges
rng = np.random.default_rng(0)
for G in [G1, G2, G3]:
    for _ in range(15):
        u, v = rng.integers(0, N, 2)
        if u != v: G.add_edge(int(u), int(v))

all_nodes = list(range(N))
result    = greedy_mcis(G1, G2, G3, all_nodes, seed=0)
is_iso    = verify_isomorphic(G1, G2, G3, result)
edges_in  = len([(u,v) for u,v in G1.edges() if u in result and v in result])

print(f"  Ground truth MCIS size: {K}")
print(f"  Found MCIS size:        {len(result)}")
print(f"  Isomorphic:             {is_iso}")
print(f"  Consensus edges found:  {edges_in}")
print(f"  MCIS nodes:             {sorted(result)}")
print(f"  Ground truth nodes:     {consensus_nodes}")
print(f"  Overlap:                {len(result & set(consensus_nodes))}/{K}")


# ── Example 2: Scaling on random Erdős–Rényi graphs ────────────────

print("\n" + "=" * 55)
print("Example 2: Scaling on random G(n, p) graphs")
print("=" * 55)

import time

sizes   = [20, 50, 100, 200]
results = []
for n in sizes:
    p = 3.0 / n  # sparse, ~3 edges per node
    G1r = nx.gnp_random_graph(n, p, directed=True, seed=1)
    G2r = nx.gnp_random_graph(n, p, directed=True, seed=2)
    G3r = nx.gnp_random_graph(n, p, directed=True, seed=3)
    t0  = time.time()
    res = greedy_mcis(G1r, G2r, G3r, list(range(n)), seed=0)
    rt  = time.time() - t0
    results.append({'n': n, 'mcis_n': len(res), 'runtime': rt})
    print(f"  n={n:4d}: MCIS={len(res):3d}, time={rt:.3f}s")


# ── Visualize ──────────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#0D1117')

# Panel 1: Example 1 circuit
ax1.set_facecolor('#0D1117')
pos = nx.circular_layout(G1)
nx.draw_networkx_nodes(G1, pos, ax=ax1,
    nodelist=list(result),
    node_color='#FFD700', node_size=250, alpha=0.95, label='MCIS')
nx.draw_networkx_nodes(G1, pos, ax=ax1,
    nodelist=[n for n in range(N) if n not in result],
    node_color='#333333', node_size=100, alpha=0.5, label='Not in MCIS')
nx.draw_networkx_edges(G1, pos, ax=ax1,
    edgelist=[(u,v) for u,v in consensus_edges if u in result and v in result],
    edge_color='#4FC3F7', arrows=True, arrowsize=15, width=2, alpha=0.9,
    connectionstyle='arc3,rad=0.15')
nx.draw_networkx_edges(G1, pos, ax=ax1,
    edgelist=[(u,v) for u,v in G1.edges() if not ((u,v) in consensus_edges)],
    edge_color='#555555', arrows=True, arrowsize=8, width=0.8, alpha=0.4)
nx.draw_networkx_labels(G1, pos, ax=ax1, font_color='white', font_size=8)
ax1.legend(facecolor='#161B22', labelcolor='white', fontsize=9)
ax1.set_title(f'Toy Example: MCIS={len(result)}/{N} nodes found\n(ground truth = {K}, ring topology)',
              color='white', fontsize=11, pad=10)
ax1.axis('off')

# Panel 2: Scaling
ax2.set_facecolor('#161B22')
ns   = [r['n'] for r in results]
rts  = [r['runtime'] for r in results]
mcis = [r['mcis_n'] for r in results]
ax2_r = ax2.twinx(); ax2_r.set_facecolor('#161B22')
l1, = ax2.plot(ns, rts, color='#4FC3F7', lw=2.5, marker='o', ms=7, label='Runtime (s)')
l2, = ax2_r.plot(ns, mcis, color='#FFD700', lw=2.5, marker='s', ms=7, ls='--', label='MCIS size')
ax2.set_xlabel('Number of nodes N', color='#aaa', fontsize=11)
ax2.set_ylabel('Runtime (seconds)', color='#4FC3F7', fontsize=11)
ax2_r.set_ylabel('MCIS size', color='#FFD700', fontsize=11)
ax2.tick_params(colors='#4FC3F7'); ax2_r.tick_params(colors='#FFD700')
ax2.spines[:].set_color('#333'); ax2_r.spines[:].set_color('#333')
ax2.set_title('Algorithm Scaling on Random Graphs', color='white', fontsize=11, pad=10)
ax2.legend(handles=[l1, l2], facecolor='#1a1a2e', labelcolor='white', fontsize=9)

fig.suptitle('mcis_connectome: Toy Demo — Algorithm Verification & Scaling',
             color='white', fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig('examples/toy_demo_output.png', dpi=150, bbox_inches='tight', facecolor='#0D1117')
print("\nSaved examples/toy_demo_output.png")
print("\nTo use with real FlyWire data:")
print("  from src.mcis_connectome import MCISSolver")
print("  solver = MCISSolver(edge_lists={...}, triplets_path='banc_meta.feather')")
print("  result = solver.solve(n_seeds=20)")
print("  print(result.summary())")
