# Examples

## `toy_mcis_demo.py` — Algorithm verification & scaling (no large data needed)

Demonstrates the MCIS algorithm on small synthetic graphs where the ground truth is known.  
Runs in seconds on any laptop. No large data downloads required.

```bash
python examples/toy_mcis_demo.py
```

**What it does:**
- Example 1: 20-node graph with a planted 8-node ring circuit → verifies the algorithm recovers it
- Example 2: Scales across random G(n,p) graphs (n = 20, 50, 100, 200) → runtime benchmark

**Output:** `examples/toy_demo_output.png`

## Using with real FlyWire data

```python
from src.mcis_connectome import MCISSolver

solver = MCISSolver(
    edge_lists={
        'BANC': 'path/to/banc_626_edge_list.csv',
        'FAFB': 'path/to/fafb_783_edge_list.csv',
        'MANC': 'path/to/manc_1.2.1_edge_list.csv',
    },
    triplets_path='path/to/banc_meta.feather',  # download URL in main README
    n_seeds=20,
)
result = solver.solve()
print(result.summary())
# MCISResult(N=109, edges=14, isomorphic=True)
result.to_csv('network.csv')
```

Or via CLI:
```bash
python -m mcis_connectome.cli \
    --banc  banc_626_edge_list.csv \
    --fafb  fafb_783_edge_list.csv \
    --manc  manc_1.2.1_edge_list.csv \
    --meta  banc_meta.feather \
    --seeds 20 \
    --out   network.csv
```
