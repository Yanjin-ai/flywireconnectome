"""
FlyWire Cross-Connectome Maximum Common Induced Subgraph Analysis
=================================================================
Goal: Find the largest neuronal circuit (directed induced subgraph)
      shared across at least 3 of 5 connectome datasets.

Strategy: Use cell type annotations as neuron correspondence keys.
          Build cell-type-level graphs, find maximum common induced subgraph.

Author: Analysis pipeline for FlyWire Qualification Challenge
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import time
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = '/Volumes/SANDISK ELE/flywire研究/'

# ============================================================
# STEP 1: Load Edge Lists
# ============================================================

def load_edge_list(filepath, name):
    """Load a connectome edge list CSV into a NetworkX DiGraph."""
    print(f"  Loading {name}...")
    t0 = time.time()
    df = pd.read_csv(filepath, dtype=str)
    df.columns = ['source', 'target']
    G = nx.from_pandas_edgelist(df, source='source', target='target',
                                 create_using=nx.DiGraph())
    print(f"    → {G.number_of_nodes():,} neurons, {G.number_of_edges():,} synaptic edges ({time.time()-t0:.1f}s)")
    return G

print("=" * 60)
print("STEP 1: Loading connectome edge lists")
print("=" * 60)

graphs = {}
graphs['FAFB'] = load_edge_list(DATA_DIR + 'fafb_783_edge_list.csv', 'FAFB')
graphs['BANC'] = load_edge_list(DATA_DIR + 'banc_626_edge_list (2).csv', 'BANC')
graphs['MANC'] = load_edge_list(DATA_DIR + 'manc_1.2.1_edge_list.csv', 'MANC')
graphs['MAOL'] = load_edge_list(DATA_DIR + 'maol_1.1_edge_list.csv', 'MAOL')
graphs['MCNS'] = load_edge_list(DATA_DIR + 'mcns_0.9_edge_list.csv', 'MCNS')

# ============================================================
# STEP 2: Load Cell Type Annotations
# ============================================================

print("\n" + "=" * 60)
print("STEP 2: Loading cell type annotations")
print("=" * 60)

# FAFB: root_id → cell_type
fafb_ann = pd.read_table(DATA_DIR + 'fafb_annotations.tsv', low_memory=False,
                          usecols=['root_id','cell_type','super_class','side','top_nt'])
fafb_ann['root_id'] = fafb_ann['root_id'].astype(str)
fafb_ann = fafb_ann.dropna(subset=['cell_type'])
fafb_id2type = dict(zip(fafb_ann['root_id'], fafb_ann['cell_type']))
print(f"  FAFB: {len(fafb_id2type):,} neurons with cell_type, {fafb_ann['cell_type'].nunique():,} unique types")

# MCNS: bodyId → flywireType (same naming convention as FAFB)
mcns_ann = pd.read_feather(DATA_DIR + 'mcns_annotations.feather',
                            columns=['bodyId','flywireType','type','superclass','mancBodyid','mancType','somaSide'])
mcns_ann['bodyId'] = mcns_ann['bodyId'].astype(str)
mcns_ann = mcns_ann.dropna(subset=['flywireType'])
mcns_id2type = dict(zip(mcns_ann['bodyId'], mcns_ann['flywireType']))
print(f"  MCNS: {len(mcns_id2type):,} neurons with flywireType, {mcns_ann['flywireType'].nunique():,} unique types")

# MANC: use mancBodyid from MCNS annotations (MCNS has cross-links to MANC)
manc_link = mcns_ann[mcns_ann['mancBodyid'].notna()][['bodyId','flywireType','mancBodyid','mancType']].copy()
manc_link['mancBodyid'] = manc_link['mancBodyid'].astype(float).astype(int).astype(str)
# Use flywireType as the canonical type for MANC neurons too
manc_id2type = dict(zip(manc_link['mancBodyid'], manc_link['flywireType']))
# Fill with mancType where flywireType is missing
for row in manc_link[manc_link['flywireType'].isna()].itertuples():
    if pd.notna(row.mancType):
        manc_id2type[row.mancBodyid] = row.mancType
print(f"  MANC: {len(manc_id2type):,} neurons with type (via MCNS cross-links)")

# MAOL: load optic lobe matching table
maol_ann = pd.read_table(DATA_DIR + 'maol_annotations.feather')  # actually TSV
maol_ann = maol_ann.rename(columns={'OL_type':'maol_type', 'Schlegel_type':'schlegel_type'})
print(f"  MAOL: loaded {len(maol_ann):,} optic lobe type matches")
print(f"  MAOL columns: {list(maol_ann.columns)}")

# ============================================================
# STEP 3: Find Shared Cell Types Across Dataset Pairs
# ============================================================

print("\n" + "=" * 60)
print("STEP 3: Computing shared cell types across datasets")
print("=" * 60)

fafb_types = set(fafb_id2type.values())
mcns_types = set(mcns_id2type.values())
manc_types = set(manc_id2type.values())

overlap_fafb_mcns = fafb_types & mcns_types
overlap_fafb_manc = fafb_types & manc_types
overlap_mcns_manc = mcns_types & manc_types
overlap_all3 = fafb_types & mcns_types & manc_types

print(f"  FAFB ∩ MCNS: {len(overlap_fafb_mcns):,} shared cell types")
print(f"  FAFB ∩ MANC: {len(overlap_fafb_manc):,} shared cell types")
print(f"  MCNS ∩ MANC: {len(overlap_mcns_manc):,} shared cell types")
print(f"  FAFB ∩ MCNS ∩ MANC: {len(overlap_all3):,} shared cell types")

# ============================================================
# STEP 4: Build Cell-Type-Level Graphs
# ============================================================

print("\n" + "=" * 60)
print("STEP 4: Building cell-type-level graphs")
print("=" * 60)

def build_celltype_graph(G_neuron, id2type, name):
    """
    Collapse a neuron-level graph into a cell-type-level graph.

    An edge A→B exists in the cell-type graph if ANY neuron of type A
    synapses onto ANY neuron of type B in the neuron graph.

    This is the key compression step: ~100k neurons → ~8k cell types
    """
    CT_G = nx.DiGraph()
    mapped = 0
    unmapped = 0

    for u, v in G_neuron.edges():
        type_u = id2type.get(u)
        type_v = id2type.get(v)
        if type_u and type_v and type_u != type_v:
            CT_G.add_edge(type_u, type_v)
            mapped += 1
        else:
            unmapped += 1

    print(f"  {name}: {CT_G.number_of_nodes():,} cell types, {CT_G.number_of_edges():,} type-level edges")
    print(f"    ({mapped:,} neuron edges mapped, {unmapped:,} unmapped/self-type)")
    return CT_G

ct_graphs = {}
ct_graphs['FAFB'] = build_celltype_graph(graphs['FAFB'], fafb_id2type, 'FAFB')
ct_graphs['MCNS'] = build_celltype_graph(graphs['MCNS'], mcns_id2type, 'MCNS')
ct_graphs['MANC'] = build_celltype_graph(graphs['MANC'], manc_id2type, 'MANC')

# ============================================================
# STEP 5: Find Maximum Common Induced Subgraph (MCIS)
# ============================================================

print("\n" + "=" * 60)
print("STEP 5: Searching for Maximum Common Induced Subgraph")
print("=" * 60)

def get_induced_subgraph(G, node_set):
    """Return the induced subgraph: only edges between nodes in node_set."""
    return G.subgraph(node_set).copy()

def check_isomorphic_induced(G1, G2, G3, node_set):
    """
    Check if the induced subgraphs on node_set are mutually isomorphic
    across three cell-type graphs.

    Since nodes are cell types (labeled), isomorphism means:
    the exact same edges exist between the exact same cell types.
    This is O(N²) comparison, not general NP-hard isomorphism!
    """
    sg1 = get_induced_subgraph(G1, node_set)
    sg2 = get_induced_subgraph(G2, node_set)
    sg3 = get_induced_subgraph(G3, node_set)

    # With labeled nodes, isomorphism = same edge set
    edges1 = set(sg1.edges())
    edges2 = set(sg2.edges())
    edges3 = set(sg3.edges())

    return edges1 == edges2 == edges3, edges1, edges2, edges3

def greedy_mcis_search(G1, G2, G3, shared_types, max_iter=1000):
    """
    Greedy algorithm to find Maximum Common Induced Subgraph.

    Strategy: Start with all shared cell types, greedily remove nodes
    that create edge disagreements between the three graphs.

    Returns the largest found common induced subgraph node set.
    """
    print(f"  Starting with {len(shared_types)} shared cell types...")

    # Only keep cell types that appear in ALL THREE graphs
    candidate_types = shared_types & set(G1.nodes()) & set(G2.nodes()) & set(G3.nodes())
    print(f"  Present in all 3 graphs: {len(candidate_types)} types")

    current_set = set(candidate_types)
    best_set = set()

    for iteration in range(max_iter):
        is_iso, e1, e2, e3 = check_isomorphic_induced(G1, G2, G3, current_set)

        if is_iso:
            if len(current_set) > len(best_set):
                best_set = set(current_set)
                print(f"  ✓ Iteration {iteration}: Found isomorphic subgraph with N={len(best_set)}")
            break

        # Find disagreement: edges in some graphs but not others
        all_edges = e1 | e2 | e3
        disagreement_edges = (e1 ^ e2) | (e2 ^ e3) | (e1 ^ e3)

        # Count how many disagreements each node participates in
        node_disagreements = defaultdict(int)
        for u, v in disagreement_edges:
            if u in current_set:
                node_disagreements[u] += 1
            if v in current_set:
                node_disagreements[v] += 1

        if not node_disagreements:
            best_set = set(current_set)
            break

        # Remove the node with most disagreements
        worst_node = max(node_disagreements, key=node_disagreements.get)
        current_set.remove(worst_node)

        if iteration % 100 == 0:
            print(f"  Iteration {iteration}: {len(current_set)} types remaining, "
                  f"{len(disagreement_edges)} disagreement edges")

    return best_set

# Run search on FAFB × MCNS × MANC (most biologically interesting trio)
print("\nSearching FAFB × MCNS × MANC...")
best_nodes = greedy_mcis_search(
    ct_graphs['FAFB'], ct_graphs['MCNS'], ct_graphs['MANC'],
    overlap_all3
)

print(f"\n  RESULT: Largest found common induced subgraph: N = {len(best_nodes)}")
print(f"  Cell types: {list(best_nodes)[:20]}...")

# ============================================================
# STEP 6: Build Network CSV Output
# ============================================================

print("\n" + "=" * 60)
print("STEP 6: Building network.csv")
print("=" * 60)

# For each matched cell type, find representative neuron IDs in each dataset
def get_representative_neuron(cell_type, id2type, graph):
    """Get the neuron ID in a dataset for a given cell type."""
    candidates = [nid for nid, ct in id2type.items() if ct == cell_type and nid in graph.nodes()]
    return candidates[0] if candidates else None

rows = []
for ct in best_nodes:
    fafb_nid = get_representative_neuron(ct, fafb_id2type, graphs['FAFB'])
    mcns_nid = get_representative_neuron(ct, mcns_id2type, graphs['MCNS'])
    manc_nid = get_representative_neuron(ct, manc_id2type, graphs['MANC'])
    if fafb_nid and mcns_nid and manc_nid:
        rows.append({'FAFB': fafb_nid, 'MCNS': mcns_nid, 'MANC': manc_nid, 'cell_type': ct})

output_df = pd.DataFrame(rows)
print(f"  Matched neurons: {len(output_df)} rows")

# Save network.csv (3 columns as required)
network_csv = output_df[['FAFB', 'MCNS', 'MANC']]
network_csv.to_csv(DATA_DIR + 'network.csv', index=False)
print(f"  Saved: network.csv")
print(network_csv.head(10).to_string())

# Also save with cell_type for analysis
output_df.to_csv(DATA_DIR + 'network_with_types.csv', index=False)

# ============================================================
# STEP 7: Analyze the Circuit
# ============================================================

print("\n" + "=" * 60)
print("STEP 7: Circuit Analysis")
print("=" * 60)

if len(best_nodes) > 0:
    # Get the actual induced subgraph
    circuit = get_induced_subgraph(ct_graphs['FAFB'], best_nodes)

    print(f"  Circuit size: {circuit.number_of_nodes()} nodes, {circuit.number_of_edges()} edges")

    # Degree statistics
    in_deg = dict(circuit.in_degree())
    out_deg = dict(circuit.out_degree())

    print(f"  Max in-degree: {max(in_deg.values()) if in_deg else 0}")
    print(f"  Max out-degree: {max(out_deg.values()) if out_deg else 0}")

    # Find hub nodes (high degree)
    hubs = sorted(circuit.nodes(), key=lambda n: in_deg.get(n,0)+out_deg.get(n,0), reverse=True)[:10]
    print(f"  Top hub cell types: {hubs}")

    # Annotate with superclass info
    type2superclass = dict(zip(fafb_ann['cell_type'], fafb_ann['super_class']))
    type2nt = dict(zip(fafb_ann['cell_type'], fafb_ann['top_nt']))

    superclass_dist = defaultdict(int)
    nt_dist = defaultdict(int)
    for n in best_nodes:
        superclass_dist[type2superclass.get(n, 'unknown')] += 1
        nt_dist[type2nt.get(n, 'unknown')] += 1

    print(f"  Superclass distribution: {dict(superclass_dist)}")
    print(f"  Neurotransmitter distribution: {dict(nt_dist)}")

    # Check for known circuit motifs
    triangles = sum(nx.triangles(circuit.to_undirected()).values()) // 3
    print(f"  Triangles (feedforward loops candidate): {triangles}")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
print(f"  network.csv saved with {len(output_df)} matched neuron pairs")
print(f"  Circuit N = {len(best_nodes)} cell types")
print(f"  Datasets used: FAFB (female brain) × MCNS (male full CNS) × MANC (male nerve cord)")
