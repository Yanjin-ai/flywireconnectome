"""
Command-line interface for mcis_connectome.

Usage:
    python -m mcis_connectome.cli \\
        --banc  banc_626_edge_list.csv \\
        --fafb  fafb_783_edge_list.csv \\
        --manc  manc_1.2.1_edge_list.csv \\
        --meta  banc_meta.feather \\
        --seeds 100 \\
        --out   network.csv
"""
import argparse
from .solver import MCISSolver


def main():
    parser = argparse.ArgumentParser(
        description='Find Maximum Common Induced Subgraph across BANC/FAFB/MANC connectomes.')
    parser.add_argument('--banc',  required=True, help='BANC edge list CSV')
    parser.add_argument('--fafb',  required=True, help='FAFB edge list CSV')
    parser.add_argument('--manc',  required=True, help='MANC edge list CSV')
    parser.add_argument('--meta',  required=True, help='BANC metadata feather file')
    parser.add_argument('--seeds', type=int, default=10, help='Number of random seeds (default: 10)')
    parser.add_argument('--out',   default='network.csv', help='Output CSV path (default: network.csv)')
    args = parser.parse_args()

    solver = MCISSolver(
        edge_lists={'BANC': args.banc, 'FAFB': args.fafb, 'MANC': args.manc},
        triplets_path=args.meta,
        n_seeds=args.seeds
    )
    result = solver.solve()
    print(result.summary())
    result.to_csv(args.out)


if __name__ == '__main__':
    main()
