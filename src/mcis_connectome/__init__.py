"""
mcis_connectome: Maximum Common Induced Subgraph for connectome data
====================================================================
Finds the largest set of neurons with identical directed connectivity
across multiple connectome datasets.

Usage:
    from mcis_connectome import MCISSolver
    solver = MCISSolver(edge_lists={'BANC': 'banc.csv', 'FAFB': 'fafb.csv', 'MANC': 'manc.csv'},
                        triplets='banc_meta.feather')
    result = solver.solve(n_seeds=100)
    result.to_csv('network.csv')
"""
from .solver import MCISSolver
from .utils import load_edge_list, build_consensus_component

__version__ = '1.0.0'
__all__ = ['MCISSolver', 'load_edge_list', 'build_consensus_component']
