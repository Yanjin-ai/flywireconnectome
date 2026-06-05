"""
RECONSTRUCTED EXPERIMENT: Cross-Connectome MCIS Analysis
=========================================================

PREVIOUS APPROACH (WRONG):
- Cell type level matching → aggregated connectivity → missed neuron-level structure
- Self-invented greedy without understanding data
- Never used BANC despite it being the KEY dataset
- No literature grounding

THIS APPROACH (CORRECT, LITERATURE-GROUNDED):
- Uses BANC meta file's fafb_match + manc_match columns
  → INDIVIDUAL NEURON LEVEL 3-way correspondence
  → Established by FlyWire team via NBLAST morphological matching
  → The SAME method as Schlegel 2024, Bates 2025
- Dataset triplet: BANC × FAFB × MANC
  → BANC edge list uses root_626 IDs
  → fafb_match gives root_id in fafb_783_edge_list.csv
  → manc_match gives bodyId in manc_1.2.1_edge_list.csv
- Algorithm: MCIS via labeled-node induced subgraph comparison
  → With neuron-level bijection, isomorphism = edge set equality
  → Greedy + expansion search (reproducible, seed-fixed)

THEORETICAL BASIS:
- Witvliet et al. 2021: connectome stereotypy methodology
- Schlegel et al. 2024: NBLAST cross-connectome cell typing
- Bates et al. 2025: BANC/FAFB/MANC comparative analysis
- Sensorimotor bottleneck hypothesis (Pospisil et al. 2024)

ASSUMPTIONS (explicitly stated):
1. One representative neuron per cell type when type has multiple neurons
   → We use all matched neurons (many types have 2: L+R)
2. Edge existence (not weight) defines the circuit structure
   → Per challenge specification
3. NBLAST-based correspondence from BANC meta is ground truth
   → Established by expert annotation, not our assumption
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import time
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = '/Volumes/SANDISK ELE/flywire研究/'
np.random.seed(42)

print("=" * 65)
print("STEP 1: Load BANC Meta — the 3-way neuron correspondence table")
print("=" * 65)

banc_meta = pd.read_feather(DATA_DIR + 'banc_meta.feather')
print(f"  BANC meta: {len(banc_meta):,} neurons, {banc_meta.shape[1]} columns")

# Build the 3-way matched set
# BANC uses root_626 (matches banc_626_edge_list.csv)
# fafb_match → root_id in fafb_783_edge_list.csv
# manc_match → bodyId in manc_1.2.1_edge_list.csv
matched = banc_meta[
    banc_meta['fafb_match'].notna() &
    banc_meta['manc_match'].notna()
][['root_626', 'fafb_match', 'manc_match',
   'cell_type', 'fafb_cell_type', 'manc_cell_type',
   'super_class', 'flow', 'neurotransmitter_predicted',
   'sexually_dimorphic']].copy()

matched['root_626'] = matched['root_626'].astype(str)
matched['fafb_match'] = matched['fafb_match'].astype(str)
matched['manc_match'] = matched['manc_match'].astype(str)
matched = matched.drop_duplicates(subset=['root_626'])

print(f"  Neurons with BOTH fafb_match and manc_match: {len(matched):,}")
print(f"  Superclass breakdown:")
print(matched['super_class'].value_counts().head(8).to_string())

# Build index: banc_id → (fafb_id, manc_id, metadata)
# Also build reverse indices
banc2fafb = dict(zip(matched['root_626'], matched['fafb_match']))
banc2manc = dict(zip(matched['root_626'], matched['manc_match']))
banc_nodes = set(matched['root_626'].tolist())

print(f"\n  Correspondence established for {len(banc_nodes):,} neurons")
print(f"  Method: NBLAST morphological matching (Schlegel 2024 methodology)")


print("\n" + "=" * 65)
print("STEP 2: Load Edge Lists — build neuron-level directed graphs")
print("=" * 65)

def load_graph_fast(filepath, node_filter=None):
    """Load edge list; optionally filter to only keep edges within node_filter."""
    t0 = time.time()
    df = pd.read_csv(filepath, dtype=str)
    df.columns = ['source', 'target']
    if node_filter is not None:
        df = df[df['source'].isin(node_filter) & df['target'].isin(node_filter)]
    G = nx.from_pandas_edgelist(df, source='source', target='target',
                                 create_using=nx.DiGraph())
    print(f"  → {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges  ({time.time()-t0:.1f}s)")
    return G

# Only load edges between matched neurons (much faster, memory efficient)
fafb_matched_ids = set(matched['fafb_match'].tolist())
manc_matched_ids = set(matched['manc_match'].tolist())
banc_matched_ids = banc_nodes

print("  BANC (edges within matched neurons):")
G_banc = load_graph_fast(DATA_DIR + 'banc_626_edge_list (2).csv', banc_matched_ids)

print("  FAFB (edges within matched neurons):")
G_fafb = load_graph_fast(DATA_DIR + 'fafb_783_edge_list.csv', fafb_matched_ids)

print("  MANC (edges within matched neurons):")
G_manc = load_graph_fast(DATA_DIR + 'manc_1.2.1_edge_list.csv', manc_matched_ids)


print("\n" + "=" * 65)
print("STEP 3: Build Unified Node Space")
print("=" * 65)

# Key insight: each BANC neuron is the "anchor"
# A matched triple is (banc_id, fafb_id, manc_id)
# The MCIS problem: find the largest subset T of triples such that
#   for any two triples (b,f,m) and (b',f',m') in T:
#   edge(b→b') in BANC ⟺ edge(f→f') in FAFB ⟺ edge(m→m') in MANC
#   AND
#   edge(b'→b) in BANC ⟺ edge(f'→f) in FAFB ⟺ edge(m'→m) in MANC

def build_adjacency_set(G, node_list, node2idx):
    """
    Build set of (i,j) index pairs representing edges,
    where i,j are indices into node_list.
    This makes cross-graph comparison O(E) not O(N²).
    """
    edges = set()
    for u, v in G.edges():
        if u in node2idx and v in node2idx:
            edges.add((node2idx[u], node2idx[v]))
    return edges

def get_induced_edges_indexed(G_nodes_set, edge_set_indexed, active_idx):
    """Get edges where both endpoints are in active_idx (as set of indices)."""
    return {(i,j) for i,j in edge_set_indexed if i in active_idx and j in active_idx}

# Build index: position in matched list → all three dataset IDs
triples = matched[['root_626','fafb_match','manc_match']].reset_index(drop=True)
n_triples = len(triples)
print(f"  Total matched triples: {n_triples:,}")

# Build indexed edge sets (i,j) where i,j are indices into triples
banc_idx = {row.root_626: idx for idx, row in triples.iterrows()}
fafb_idx = {row.fafb_match: idx for idx, row in triples.iterrows()}
manc_idx = {row.manc_match: idx for idx, row in triples.iterrows()}

print("  Building indexed edge sets...")
banc_edges_idx = set()
for u, v in G_banc.edges():
    if u in banc_idx and v in banc_idx:
        banc_edges_idx.add((banc_idx[u], banc_idx[v]))

fafb_edges_idx = set()
for u, v in G_fafb.edges():
    if u in fafb_idx and v in fafb_idx:
        fafb_edges_idx.add((fafb_idx[u], fafb_idx[v]))

manc_edges_idx = set()
for u, v in G_manc.edges():
    if u in manc_idx and v in manc_idx:
        manc_edges_idx.add((manc_idx[u], manc_idx[v]))

print(f"  BANC internal edges (matched neurons only): {len(banc_edges_idx):,}")
print(f"  FAFB internal edges (matched neurons only): {len(fafb_edges_idx):,}")
print(f"  MANC internal edges (matched neurons only): {len(manc_edges_idx):,}")

# How many edges are agreed on across all three?
all3_agree = banc_edges_idx & fafb_edges_idx & manc_edges_idx
any_dataset = banc_edges_idx | fafb_edges_idx | manc_edges_idx
disagreements = any_dataset - all3_agree
print(f"  Edges present in ALL 3 datasets: {len(all3_agree):,}")
print(f"  Edges present in ANY dataset: {len(any_dataset):,}")
print(f"  Disagreement edges: {len(disagreements):,}")
print(f"  Consensus rate: {len(all3_agree)/max(len(any_dataset),1)*100:.1f}%")


print("\n" + "=" * 65)
print("STEP 4: MCIS Search — Greedy Disagreement Removal + Expansion")
print("=" * 65)
print()
print("  ALGORITHM DESCRIPTION:")
print("  ─────────────────────")
print("  Problem: find largest set S of neuron indices such that")
print("    ∀ i,j ∈ S: (i→j) ∈ BANC ⟺ (i→j) ∈ FAFB ⟺ (i→j) ∈ MANC")
print()
print("  Strategy: Greedy Disagreement Removal")
print("    1. Start with all N triples")
print("    2. Find 'disagreement edges': in some datasets but not others")
print("    3. Remove neuron involved in most disagreements")
print("    4. Repeat until zero disagreements")
print("    5. Expansion: try adding removed neurons back (local search)")
print()
print("  Complexity: O(N × D) per iteration, N=nodes, D=disagreement edges")
print("  Reproducibility: fixed seed 42, deterministic")
print()

def mcis_greedy_remove(all_edge_sets, n_nodes, max_iter=5000):
    """
    Greedy disagreement removal to find MCIS.

    all_edge_sets: list of sets of (i,j) indexed edges
    n_nodes: total number of nodes
    Returns: set of node indices forming common induced subgraph
    """
    active = set(range(n_nodes))
    best = set()

    for iteration in range(max_iter):
        # Get induced edge sets for active nodes
        induced = [
            {(i,j) for i,j in es if i in active and j in active}
            for es in all_edge_sets
        ]

        # Check consensus
        union = set().union(*induced)
        disagree = set()
        for e in union:
            vals = [e in ind for ind in induced]
            if not all(vals):
                disagree.add(e)

        if not disagree:
            if len(active) > len(best):
                best = set(active)
                print(f"  ✓ iter {iteration:4d}: ISOMORPHIC, N={len(best)}")
            break

        # Score nodes by disagreement count
        node_score = defaultdict(int)
        for u, v in disagree:
            if u in active: node_score[u] += 1
            if v in active: node_score[v] += 1

        # Remove worst nodes (top-5 per iteration for speed)
        worst = sorted(node_score, key=node_score.get, reverse=True)[:5]
        for w in worst:
            if w in active:
                active.remove(w)

        if iteration % 500 == 0 and iteration > 0:
            print(f"  iter {iteration:4d}: {len(active):,} nodes, {len(disagree):,} disagreements")

    return best

def mcis_expand(best_set, all_edge_sets, candidate_pool):
    """
    Expansion phase: try adding each removed node back.
    A node can rejoin if it doesn't create new disagreements.
    """
    current = set(best_set)
    candidates = set(candidate_pool) - current
    improved = True
    n_rounds = 0

    while improved:
        improved = False
        n_rounds += 1
        added_this_round = []

        for node in list(candidates):
            test = current | {node}
            induced = [
                {(i,j) for i,j in es if i in test and j in test}
                for es in all_edge_sets
            ]
            union = set().union(*induced)
            has_disagree = any(
                not all(e in ind for ind in induced)
                for e in union
            )
            if not has_disagree:
                current.add(node)
                candidates.remove(node)
                added_this_round.append(node)
                improved = True

        if added_this_round:
            print(f"  Expansion round {n_rounds}: +{len(added_this_round)} nodes → N={len(current)}")

    return current

edge_sets = [banc_edges_idx, fafb_edges_idx, manc_edges_idx]

# Multi-start search for better results
print("  === MULTI-START GREEDY SEARCH ===")
global_best = set()

for trial in range(5):
    if trial == 0:
        start = set(range(n_triples))
    else:
        # Random subsample to escape local optima
        remove_n = n_triples // 4
        rng = np.random.default_rng(trial * 17)
        remove_idx = set(rng.choice(n_triples, remove_n, replace=False).tolist())
        start = set(range(n_triples)) - remove_idx

    result = mcis_greedy_remove(edge_sets, n_triples, max_iter=3000)
    # Override active set for non-zero trials
    if trial > 0:
        # Re-run from perturbed start
        active = start.copy()
        for _ in range(3000):
            induced = [{(i,j) for i,j in es if i in active and j in active} for es in edge_sets]
            union = set().union(*induced)
            disagree = {e for e in union if not all(e in ind for ind in induced)}
            if not disagree:
                result = set(active)
                break
            node_score = defaultdict(int)
            for u,v in disagree:
                if u in active: node_score[u] += 1
                if v in active: node_score[v] += 1
            for w in sorted(node_score, key=node_score.get, reverse=True)[:5]:
                if w in active: active.remove(w)
        else:
            result = set(active)

    if len(result) > len(global_best):
        global_best = result
        print(f"  Trial {trial}: NEW BEST N={len(global_best)}")
    else:
        print(f"  Trial {trial}: N={len(result)}")

print(f"\n  Post-greedy best: N={len(global_best)}")
print("  === EXPANSION PHASE ===")
final_nodes = mcis_expand(global_best, edge_sets, range(n_triples))
print(f"\n  Final N={len(final_nodes)}")

# Verify
induced_final = [{(i,j) for i,j in es if i in final_nodes and j in final_nodes} for es in edge_sets]
union_final = set().union(*induced_final)
disagree_final = {e for e in union_final if not all(e in ind for ind in induced_final)}
print(f"  Verification: disagreements = {len(disagree_final)} (should be 0)")
print(f"  Induced edges in circuit: {len(induced_final[0])}")

print("\n" + "=" * 65)
print("STEP 5: Analyze & Output")
print("=" * 65)

# Map back to neuron IDs and metadata
result_triples = triples.iloc[sorted(final_nodes)].copy()
result_triples['circuit_node_idx'] = sorted(final_nodes)

print(f"\n  Circuit size: {len(result_triples)} neurons")
print(f"  Internal edges: {len(induced_final[0])}")
print()
print("  Superclass composition:")
sc = result_triples['super_class'].value_counts()
print(sc.to_string())
print()
print("  Neurotransmitter distribution:")
nt = result_triples['neurotransmitter_predicted'].value_counts()
print(nt.to_string())
print()
print("  Sexually dimorphic?")
print(result_triples['sexually_dimorphic'].value_counts().to_string())

# Show edge structure
print()
print("  Internal circuit edges (BANC→BANC with cross-dataset verification):")
for (i,j) in induced_final[0]:
    bi = triples.iloc[i]['root_626']
    fi = triples.iloc[i]['fafb_match']
    bj = triples.iloc[j]['root_626']
    fj = triples.iloc[j]['fafb_match']
    ci = matched[matched['root_626']==bi]['cell_type'].values[0] if len(matched[matched['root_626']==bi])>0 else '?'
    cj = matched[matched['root_626']==bj]['cell_type'].values[0] if len(matched[matched['root_626']==bj])>0 else '?'
    print(f"    {ci} ({bi[:8]}...) → {cj} ({bj[:8]}...)")

# Save network.csv — FORMAT: 3 columns, N rows, neuron IDs per dataset
output = result_triples[['root_626', 'fafb_match', 'manc_match']].copy()
output.columns = ['BANC', 'FAFB', 'MANC']
output.to_csv(DATA_DIR + 'network.csv', index=False)

# Save enriched version for analysis
result_triples.to_csv(DATA_DIR + 'network_enriched.csv', index=False)

print(f"\n  Saved: network.csv ({len(output)} rows × 3 columns)")
print(f"  Format: BANC (root_626) | FAFB (fafb_783) | MANC (manc_1.2.1)")
print()
print("  Preview:")
print(output.head(10).to_string())

print("\n" + "=" * 65)
print("STEP 6: Second search — Focus on DN/AN subpopulation")
print("=" * 65)
# Hypothesis: sensorimotor neurons (DN+AN) are most conserved
# Test: restrict to these and find larger MCIS

sm_mask = matched['super_class'].isin(['descending','ascending','sensory_ascending'])
sm_triples = matched[sm_mask].reset_index(drop=True)
n_sm = len(sm_triples)
print(f"  Sensorimotor neurons (DN+AN): {n_sm}")

if n_sm > 0:
    sm_banc_idx = {row.root_626: idx for idx, row in sm_triples.iterrows()}
    sm_fafb_idx = {row.fafb_match: idx for idx, row in sm_triples.iterrows()}
    sm_manc_idx = {row.manc_match: idx for idx, row in sm_triples.iterrows()}

    sm_banc_e = {(sm_banc_idx[u], sm_banc_idx[v]) for u,v in G_banc.edges()
                  if u in sm_banc_idx and v in sm_banc_idx}
    sm_fafb_e = {(sm_fafb_idx[u], sm_fafb_idx[v]) for u,v in G_fafb.edges()
                  if u in sm_fafb_idx and v in sm_fafb_idx}
    sm_manc_e = {(sm_manc_idx[u], sm_manc_idx[v]) for u,v in G_manc.edges()
                  if u in sm_manc_idx and v in sm_manc_idx}

    sm_edge_sets = [sm_banc_e, sm_fafb_e, sm_manc_e]

    print(f"  BANC SM internal edges: {len(sm_banc_e)}")
    print(f"  FAFB SM internal edges: {len(sm_fafb_e)}")
    print(f"  MANC SM internal edges: {len(sm_manc_e)}")

    sm_result = mcis_greedy_remove(sm_edge_sets, n_sm, max_iter=3000)
    sm_expanded = mcis_expand(sm_result, sm_edge_sets, range(n_sm))

    # Verify
    sm_induced = [{(i,j) for i,j in es if i in sm_expanded and j in sm_expanded} for es in sm_edge_sets]
    sm_union = set().union(*sm_induced)
    sm_disagree = {e for e in sm_union if not all(e in ind for ind in sm_induced)}
    print(f"  SM-specific MCIS: N={len(sm_expanded)}, edges={len(sm_induced[0])}, disagree={len(sm_disagree)}")

    if len(sm_expanded) > len(final_nodes):
        print("  ★ SM-focused search found LARGER circuit! Updating network.csv")
        sm_result_triples = sm_triples.iloc[sorted(sm_expanded)]
        sm_output = sm_result_triples[['root_626','fafb_match','manc_match']].copy()
        sm_output.columns = ['BANC','FAFB','MANC']
        sm_output.to_csv(DATA_DIR + 'network.csv', index=False)
        print(f"  Saved updated network.csv with N={len(sm_expanded)}")

print("\n" + "=" * 65)
print("COMPLETE")
print("=" * 65)
print(f"  Method: NBLAST-matched neuron triplets + MCIS greedy search")
print(f"  Datasets: BANC (root_626) × FAFB (v783) × MANC (v1.2.1)")
print(f"  N = {len(final_nodes)} verified isomorphic neurons")
print(f"  Theoretical basis: Schlegel 2024, Bates 2025, sensorimotor bottleneck")
print(f"  Reproducibility: deterministic with seed=42")
