"""Utility functions for edge list loading and graph construction."""
import pandas as pd
import networkx as nx
from typing import Set


def load_edge_list(path: str, node_filter: Set[str] = None) -> nx.DiGraph:
    """
    Load a connectome edge list CSV into a directed graph.

    Parameters
    ----------
    path : str
        Path to CSV with 'source neuron id' and 'target neuron id' columns.
    node_filter : set, optional
        If given, keep only edges where both endpoints are in this set.

    Returns
    -------
    nx.DiGraph
    """
    df = pd.read_csv(path, dtype=str)
    df.columns = ['source', 'target']
    if node_filter is not None:
        df = df[df['source'].isin(node_filter) & df['target'].isin(node_filter)]
    return nx.from_pandas_edgelist(df, 'source', 'target', create_using=nx.DiGraph())


def build_consensus_component(be, fe, me):
    """
    Build the giant weakly-connected component of the consensus (all-3-agree) graph.

    Parameters
    ----------
    be, fe, me : set of (int,int) tuples
        Indexed edge sets for each dataset.

    Returns
    -------
    giant : list of int (node indices in giant component)
    gbe, gfe, gme : edge sets restricted to giant component (local indices)
    """
    agree = be & fe & me
    CG = nx.DiGraph()
    CG.add_edges_from(agree)
    if not CG.nodes:
        return [], set(), set(), set()
    giant  = sorted(max(nx.weakly_connected_components(CG), key=len))
    loc    = {g: l for l, g in enumerate(giant)}
    gbe    = {(loc[u], loc[v]) for u,v in be if u in loc and v in loc}
    gfe    = {(loc[u], loc[v]) for u,v in fe if u in loc and v in loc}
    gme    = {(loc[u], loc[v]) for u,v in me if u in loc and v in loc}
    return giant, gbe, gfe, gme
