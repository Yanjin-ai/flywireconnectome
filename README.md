# Structural Invariance at the Sensorimotor Interface

> *Maximum Common Induced Subgraph across three independent Drosophila connectomes*  
> FlyWire Qualification Challenge · June 2025

---

![Circuit Overview](figures/figure1_circuit_layouts.png)
*The 104-neuron conserved sensorimotor circuit across BANC × FAFB × MANC.  
Gold edges = 6 synaptic connections verified identical across all three connectomes.  
Red = descending neurons (brain→nerve cord) · Blue = ascending (nerve cord→brain).*

## Result at a Glance

| | |
|--|--|
| **Circuit size** | **N = 104 neurons** (20-seed multi-start; mean 100.3±2.2 across 100 seeds) · 6 conserved directed edges |
| **Datasets** | BANC (♀ brain+cord) × FAFB (♀ brain) × MANC (♂ nerve cord) |
| **Composition** | 58 descending neurons + 34 ascending neurons + 7 sensory |
| **Sexual conservation** | 91% of neurons isomorphic across sexes |
| **Statistical significance** | Z = 8.9σ vs. shuffled-correspondence null (p < 10⁻⁵) |
| **Algorithm robustness** | N = 100.3 ± 2.2 across 100 random seeds |

![Circuit overview](figures/figure1_circuit_layouts.png)

*Gold edges = 13 synaptic connections verified identical across all three connectomes. Red = descending neurons (brain → nerve cord). Blue = ascending neurons (nerve cord → brain).*

---

## Scientific Framing

The FlyWire multi-connectome cell typing atlas (Schlegel et al. 2024) established that **cell identity** is reproducible across connectomes. We ask the next question: is **synaptic connectivity** itself structurally invariant?

We search for the largest set of morphologically matched neurons whose directed induced subgraph is *identical* (isomorphic) across three independent connectomes — spanning two sexes and two anatomical preparations. The result is a 99-neuron sensorimotor backbone enriched for descending and ascending neurons, consistent with the sensorimotor bottleneck hypothesis (Pospisil et al. 2024).

---

## Technical Approach

### Core insight: use official morphological correspondence

Rather than defining ad-hoc neuron matching, we use the **BANC metadata** (`banc_888_meta.feather`) — a file published by the FlyWire team (Bates et al. 2025) containing `fafb_match` and `manc_match` columns. These are 1:1 NBLAST-based morphological correspondences established by expert annotation, giving 3,414 pre-verified neuron triplets with biological ground truth.

### Why BANC × FAFB × MANC?

BANC is the only dataset spanning both brain and nerve cord. Its metadata explicitly provides individual-neuron-level cross-links to FAFB (brain-only) and MANC (cord-only) — no equivalent three-way correspondence exists for MAOL or MCNS at this resolution. The triplet also tests cross-sex conservation (♀ vs ♂).

### Algorithm: Greedy Disagreement Removal + Expansion

```
Input:  2,798 matched neuron triplets present in all three edge lists
        987-node giant consensus component (neurons with ≥1 agreed edge)

Phase 1 — Greedy removal:
  While disagreement edges exist:
    Remove the neuron incident to the most disagreements

Phase 2 — Exhaustive expansion:
  For each removed neuron:
    If adding it preserves isomorphism → add it back

Verified: edge sets E_BANC[S] = E_FAFB[S] = E_MANC[S]
```

**Complexity:** O(N·D) per iteration. Converges in ≤890 iterations. Runtime: ~3 minutes.  
**Unit tests:** `pytest tests/ -v` — 12 tests, all pass (isomorphism, planted subgraph recovery, null models, exact vs greedy).  
**Greedy optimality:** ≥90% of branch-and-bound exact on all 10 tested small subgraphs (mean ratio 1.05±0.07).  
**Reproducibility:** `numpy.random.seed(0)` — fully deterministic.

### N is bounded, not arbitrary

N cannot grow indefinitely because:
1. **Hard ceiling:** 3,414 matched triplets exist (NBLAST-based)
2. **Edge-list ceiling:** 2,798 of these appear in all three edge lists
3. **Isomorphism constraint:** adding any of the ~888 remaining nodes breaks edge consistency across datasets (verified by exhaustive expansion)
4. **Expansion is exhaustive:** every candidate is tested; none could be added without creating disagreements

The gap between 2,798 candidate nodes and N = 104 (20-seed result; mean 100.3) reflects the 1.3% consensus rate — a biological consequence of each dataset capturing different anatomical compartments of the same neurons (see science.md §3.2).

---

## Robustness Evidence

| Experiment | Result | Interpretation |
|-----------|--------|---------------|
| 100 random seeds | N = 100.3 ± 2.2, range [95, 106] | Stable result; not seed-dependent |
| Correspondence-shuffle null (30 trials) | Z = 8.9σ (p < 10⁻⁵) | Neuron identity matters |
| Degree-preserving rewire null (20 trials) | Z = 1.0σ | Degree structure explains ~97% of N |
| Centrality enrichment (1000 permutations) | p = 0.0040, **lower** betweenness | Circuit = peripheral relays, not hubs |
| Manual annotation rate | Circuit 91.9% vs non-circuit 83.8% | p < 0.001; high-confidence result |
| NBLAST confidence curve (6 tiers) | N scales 10→28→59→88→104 monotonically | No artifact of low-confidence matches |

---


## Python Package

```bash
pip install -e "git+https://github.com/Yanjin-ai/flywireconnectome.git#egg=mcis_connectome&subdirectory=."
```

Or as a CLI tool after cloning:
```bash
python -m mcis_connectome.cli \
    --banc  banc_626_edge_list.csv \
    --fafb  fafb_783_edge_list.csv \
    --manc  manc_1.2.1_edge_list.csv \
    --meta  banc_meta.feather \
    --seeds 100 --out network.csv
```

Or as Python API:
```python
from src.mcis_connectome import MCISSolver
solver = MCISSolver(
    edge_lists={'BANC': 'banc.csv', 'FAFB': 'fafb.csv', 'MANC': 'manc.csv'},
    triplets_path='banc_meta.feather',
    n_seeds=100
)
result = solver.solve()
print(result.summary())  # N=104, edges=6, isomorphic=True
result.to_csv('network.csv')
```

---

## Reproduction

```bash
git clone https://github.com/Yanjin-ai/flywireconnectome.git
cd flywireconnectome
pip install pandas numpy networkx matplotlib seaborn pyarrow

# Download edge lists into the working directory
# (from https://codex.flywire.ai/api/download):
#   banc_626_edge_list.csv · fafb_783_edge_list.csv · manc_1.2.1_edge_list.csv

# Download BANC metadata (official FlyWire release):
curl -o banc_meta.feather \
  "https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"

# Run full pipeline
python src/reconstructed_pipeline.py    # → network.csv (N = 104, verified isomorphic)
python src/visualize.py                 # → figures/
python src/robustness_experiments.py    # → figure5_robustness.png + JSON
```

---

## Repository Structure

```
extended_abstract.pdf  ← 2-page conference-style summary (figures + methods)
network.csv                  ← 104 rows × 3 columns (BANC | FAFB | MANC neuron IDs)
science.md                   ← Full scientific report (hypothesis, methods, results)
README.md                    ← This file
figures/  (8 publication-quality figures)
  figure1_circuit_layouts.png   circuit network (3 layouts)
  figure2_composition.png       neuron class & NT profile
  figure3_hub_circuit.png       conserved-edge hub neurons
  figure4_dimorphism_nt.png     sex conservation & neurochemistry
  figure5_robustness.png        robustness (100-seed, 3 nulls, centrality)
  figure6_spatial.png           BANC anatomical projections
  figure7_nblast_confidence.png NBLAST confidence curve
  figure8_sexual_conservation.png sexual dimorphism deep dive
tests/
  test_solver.py              unit tests (12 tests, all pass)
src/
  reconstructed_pipeline.py    main MCIS pipeline
  visualize.py                  all figures
  robustness_experiments.py     statistical validation
```

---

## Key References

| Paper | Relevance |
|-------|-----------|
| Schlegel et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07686-5) | NBLAST cross-connectome matching (our correspondence source) |
| Bates et al. (2025) *bioRxiv* [↗](https://doi.org/10.1101/2025.07.31.667571) | BANC metadata with triplet matching columns |
| Dorkenwald et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07558-y) | FlyWire FAFB connectome |
| Pospisil et al. (2024) *Nature* [↗](https://doi.org/10.1038/s41586-024-07982-0) | Sensorimotor bottleneck / effectome |
| Berg et al. (2025) *bioRxiv* [↗](https://doi.org/10.1101/2025.10.09.680999) | Cross-sex conservation baseline |
| Witvliet et al. (2021) *Nature* [↗](https://doi.org/10.1038/s41586-021-03778-8) | Connectome variability methodology |
