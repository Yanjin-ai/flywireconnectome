# Structural Invariance at the Sensorimotor Interface

> *Maximum Common Induced Subgraph across three independent Drosophila connectomes*  
> FlyWire Qualification Challenge · June 2026

---

![Circuit Overview](figures/figure1_circuit_layouts.png)
*The 104-neuron conserved sensorimotor circuit across BANC × FAFB × MANC.  
Gold edges = 6 synaptic connections verified identical across all three connectomes.  
Red = descending neurons (brain → nerve cord) · Blue = ascending neurons (nerve cord → brain).*

---

## Result at a Glance

| | |
|--|--|
| **Circuit size** | **N = 104 neurons** (20-seed multi-start; mean 100.3 ± 2.2 across 100 seeds) |
| **Conserved edges** | 6 directed edges, verified isomorphic across all 3 datasets |
| **Datasets** | BANC v626 (♀ brain+cord) × FAFB v783 (♀ brain) × MANC v1.2.1 (♂ nerve cord) |
| **Composition** | 61 descending (58.7%) + 36 ascending (34.6%) + 7 sensory neurons |
| **Sexual conservation** | 88.5% isomorphic across ♀ and ♂ (92/104) |
| **Statistical significance** | Z = 8.9σ vs correspondence-shuffle null (p < 10⁻⁵) |
| **Cell-type enrichment** | Descending 62.7× (p=1.6×10⁻⁹⁴), Ascending 27.5× (p=2.6×10⁻⁴¹) vs FAFB background |
| **Anatomical position** | Cervical connective (Fig. 4): expected locus for brain–cord relay neurons |

---

## Scientific Framing

The FlyWire multi-connectome cell typing atlas (Schlegel et al. 2024) established that **cell-type identity** is reproducible across connectomes at the morphological level. We ask the next question: is **synaptic connectivity itself** structurally invariant?

We search for the largest set of morphologically matched neurons whose directed induced subgraph is *identical* (isomorphic) across three independent connectomes — spanning two sexes and two anatomical preparations. The result is a 104-neuron sensorimotor backbone enriched 62.7× for descending neurons and 27.5× for ascending neurons relative to the whole-brain background, consistent with the sensorimotor bottleneck hypothesis (Pospisil et al. 2024).

---

## Technical Approach

### Core insight: official morphological correspondence

Rather than defining ad hoc neuron matching, we use the **BANC metadata** (`banc_888_meta.feather`, Bates et al. 2025) which contains `fafb_match` and `manc_match` columns — 1:1 NBLAST-based morphological correspondences per neuron, giving 3,414 pre-verified triplets with biological ground truth.

### Why BANC × FAFB × MANC?

BANC is the only dataset spanning both brain and ventral nerve cord. Its metadata explicitly provides individual-neuron-level cross-links to FAFB (brain-only) and MANC (cord-only). No equivalent three-way correspondence exists for MAOL or MCNS at this resolution. The triplet also tests cross-sex conservation (♀ FAFB/BANC vs ♂ MANC).

### Algorithm: Greedy Disagreement Removal + Exhaustive Expansion

```
Input:  2,798 triplets present in all 3 edge lists
        987-node giant consensus component

Phase 1 — Greedy removal:
  While disagreement edges exist:
    Remove the neuron incident to the most disagreements
    (randomised tie-breaking; fixed seeds for reproducibility)

Phase 2 — Exhaustive expansion:
  For each removed neuron:
    If adding it preserves isomorphism → add it back

Verified: E_BANC[S] = E_FAFB[S] = E_MANC[S]
```

**Complexity:** O(N·D) per iteration. Converges in ≤890 iterations (~9 seconds).  
**Reproducibility:** Fixed random seeds; fully deterministic per seed; results reported as 20-seed multi-start.  
**Unit tests:** `pytest tests/ -v` — 12 tests, all pass.  
**Near-optimality:** ≥90% of time-limited branch-and-bound on all 10 tested small subgraphs.

### Why N is bounded

N cannot grow indefinitely:
1. **NBLAST ceiling:** 3,414 matched triplets (FlyWire team's morphological matching)
2. **Edge-list ceiling:** 2,798 appear in all three edge lists
3. **Consensus component:** 987 nodes have ≥1 agreed edge across all 3 datasets
4. **Isomorphism constraint:** exhaustive expansion verifies no further neuron can be added without breaking edge-set equality

---

## Robustness Evidence

| Experiment | Result | Interpretation |
|-----------|--------|---------------|
| 100 random seeds | N = 100.3 ± 2.2, range [95, 106] | Stable; not seed-dependent |
| Correspondence-shuffle null (30 trials) | N_null = 84.7 ± 2.4; Z = 8.9σ (p < 10⁻⁵) | Neuron identity is essential |
| Degree-preserving rewire null (20 trials) | N_null = 96.9 ± 2.1; Z = 1.0σ | Degree structure explains most of achievable N |
| Centrality permutation (1000 trials) | p = 0.004; circuit has *lower* betweenness | Peripheral relays, not hubs |
| Manual annotation rate | Circuit 91.9% vs non-circuit 83.8% (p < 0.001) | High-confidence correspondences |
| NBLAST confidence curve (6 tiers) | N = 10 → 28 → 30 → 59 → 88 → 104 (monotonic) | Not driven by low-confidence matches |

---

## Python Package

```bash
pip install -e "git+https://github.com/Yanjin-ai/flywireconnectome.git#egg=mcis_connectome&subdirectory=."
```

CLI:
```bash
python -m mcis_connectome.cli \
    --banc  banc_626_edge_list.csv \
    --fafb  fafb_783_edge_list.csv \
    --manc  manc_1.2.1_edge_list.csv \
    --meta  banc_meta.feather \
    --seeds 20 --out network.csv
```

Python API:
```python
from src.mcis_connectome import MCISSolver
solver = MCISSolver(
    edge_lists={'BANC': 'banc.csv', 'FAFB': 'fafb.csv', 'MANC': 'manc.csv'},
    triplets_path='banc_meta.feather',
    n_seeds=20
)
result = solver.solve()
print(result.summary())  # MCISResult(N=104, edges=6, isomorphic=True)
result.to_csv('network.csv')
```

---

## Reproduction

```bash
git clone https://github.com/Yanjin-ai/flywireconnectome.git
cd flywireconnectome
pip install pandas numpy networkx matplotlib seaborn pyarrow scipy

# Edge lists (from https://codex.flywire.ai/api/download):
#   banc_626_edge_list.csv · fafb_783_edge_list.csv · manc_1.2.1_edge_list.csv

# BANC metadata:
curl -o banc_meta.feather \
  "https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"

# FAFB annotations (for enrichment analysis):
# https://github.com/flyconnectome/flywire_annotations

# Run:
python src/reconstructed_pipeline.py    # → network.csv (104 rows, verified isomorphic)
python src/visualize.py                 # → figures/
python src/robustness_experiments.py    # → figure5_robustness.png
pytest tests/ -v                        # → 12/12 tests pass
```

---

## Repository Structure

```
network.csv                    104 rows × 3 columns (BANC | FAFB | MANC neuron IDs)
science.md                     Scientific report (10 sections, 13 references)
extended_abstract.pdf          2-page conference-style summary
README.md                      This file
figures/
  figure1_circuit_layouts.png  Network layouts (3 force-directed)
  figure2_composition.png      Neuron class and NT composition
  figure3_hub_circuit.png      Hub neurons with conserved edges
  figure4_dimorphism_nt.png    Sexual dimorphism overview
  figure5_robustness.png       Robustness panel (7 sub-figures)
  figure6_spatial.png          BANC anatomical projections
  figure7_nblast_confidence.png NBLAST confidence curve
  figure8_sexual_conservation.png Sexual conservation deep dive
  figure9_enrichment.png       Cell-type enrichment vs FAFB background
tests/
  test_solver.py               12 unit tests (all pass)
src/
  reconstructed_pipeline.py   Main MCIS pipeline
  visualize.py                 Figure generation
  robustness_experiments.py   Statistical validation
  mcis_connectome/             Installable Python package
    solver.py                  MCISSolver + MCISResult classes
    utils.py                   load_edge_list, build_consensus_component
    cli.py                     Command-line interface
examples/
  toy_mcis_demo.py             Algorithm demo on synthetic graphs (no data needed)
setup.py                       pip install -e . support
```

---

## Key References

| Paper | Relevance |
|-------|-----------|
| Schlegel et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07686-5) | Consensus cell-type atlas; NBLAST cross-connectome matching |
| Bates et al. (2025) *bioRxiv* [↗](https://doi.org/10.1101/2025.07.31.667571) | BANC dataset and triplet correspondence metadata |
| Dorkenwald et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07558-y) | FlyWire FAFB connectome |
| Pospisil et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07982-0) | Sensorimotor bottleneck / effectome |
| Berg et al. (2025) *bioRxiv* [↗](https://doi.org/10.1101/2025.10.09.680999) | Cross-sex conservation baseline |
| Witvliet et al. (2021) *Nature* [↗](https://doi.org/10.1038/s41586-021-03778-8) | Connectome stereotypy methodology |
