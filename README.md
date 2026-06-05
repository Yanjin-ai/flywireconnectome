# FlyWire Cross-Connectome Circuit Discovery

**FlyWire Qualification Challenge** — Maximum Common Induced Subgraph across Drosophila connectomes

## Result

**N = 75** neurons verified isomorphic across **BANC × FAFB × MANC**  
**6 conserved directed edges** — identical across all three datasets  
**Circuit identity:** Sensorimotor bottleneck (descending + ascending neurons)

```
network.csv   — 75 rows × 3 columns (BANC | FAFB | MANC neuron IDs)
science.md    — Scientific summary with visualizations and hypothesis
src/          — Analysis pipeline code
figures/      — All circuit visualizations
```

---

## Technical Approach

### Core Insight

Rather than inventing an arbitrary node correspondence, we use the **BANC metadata file** (`banc_888_meta.feather`) which contains `fafb_match` and `manc_match` columns — individual neuron-level correspondences established by the FlyWire team via **NBLAST morphological similarity** (Schlegel et al. 2024). This gives us 3,414 pre-verified neuron triplets (BANC ↔ FAFB ↔ MANC) with morphological ground truth.

### Dataset Selection

| Dataset | Organism | Region | Neurons | ID type |
|---------|---------|--------|---------|---------|
| BANC (v626) | ♀ adult | Brain + nerve cord | 112,885 | 18-digit CAVE ID |
| FAFB (v783) | ♀ adult | Brain only | 138,584 | 18-digit CAVE ID |
| MANC (v1.2.1) | ♂ adult | Nerve cord only | 23,641 | Integer body ID |

**Why this triplet?**  
BANC is the only dataset that spans both brain and nerve cord. Its metadata explicitly cross-links neurons to FAFB (brain-only) and MANC (cord-only), giving the richest correspondence table. This triplet also tests **cross-sex conservation** (♀ FAFB/BANC vs ♂ MANC).

### Algorithm: MCIS via Greedy Disagreement Removal

```
Given: 3,414 matched neuron triplets (b_i, f_i, m_i)
Goal:  Find largest subset S such that
       ∀ i,j ∈ S: edge(b_i→b_j) ∈ BANC  ⟺  edge(f_i→f_j) ∈ FAFB  ⟺  edge(m_i→m_j) ∈ MANC

Algorithm:
  1. Filter to neurons present in all 3 edge lists (N = 2,798)
  2. Restrict to giant consensus component (N = 987 nodes with ≥1 agreed edge)
  3. Greedy removal: iteratively remove the node involved in most disagreements
  4. Expansion: try adding each removed node back without breaking isomorphism
  5. Final verification: assert edge set equality across all 3 induced subgraphs
```

**Reproducibility:** `numpy.random.seed(42)` — fully deterministic.  
**Complexity:** O(N × D) per iteration, N = nodes, D = disagreement edges.  
**Guarantee:** Local optimum; true global maximum is NP-hard to certify, but expansion phase mitigates.

### Key Finding: Why Only 1.3% Edge Consensus?

Of ~245,000 unique edges among matched neurons, only 2,648 (1.3%) appear in all three datasets. This is **biologically meaningful, not a failure**:

- Descending neurons' dendrites are in the brain → captured by FAFB
- Their axonal outputs are in the nerve cord → captured by MANC  
- BANC captures both, but at different synapse densities

The 2,648 consensus edges represent connections that exist at **both ends of the sensorimotor axis** — the most fundamental, obligate wiring.

---

## Reproduction

```bash
# 1. Clone and install dependencies
git clone https://github.com/Yanjin-ai/flywireconnectome.git
cd flywireconnectome
pip install pandas numpy networkx matplotlib seaborn pyarrow

# 2. Download edge lists into data/ directory
#    BANC: banc_626_edge_list.csv
#    FAFB: fafb_783_edge_list.csv
#    MANC: manc_1.2.1_edge_list.csv
#    (available at codex.flywire.ai/api/download)

# 3. Download BANC metadata (provides cross-dataset neuron matching)
curl -o data/banc_meta.feather \
  "https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"

# 4. Run pipeline
python src/reconstructed_pipeline.py   # → network.csv
python src/visualize.py                # → figures/
```

Expected output: `network.csv` with 75 rows, 3 columns.  
Runtime: ~3 minutes on a standard laptop.

---

## Assumptions

1. **Neuron correspondence** is defined by NBLAST morphological matching (Schlegel et al. 2024), not our invention — we inherit the FlyWire team's ground-truth matching.
2. **Edge definition**: an edge exists if any synapse is recorded between the two neurons (unweighted, per challenge specification).
3. **Representative neuron**: for cell types with multiple neurons (bilateral pairs), any matched representative is valid for the structural isomorphism check.
4. **Local optimality**: the greedy algorithm finds a local maximum; the true MCIS may be larger.

---

## Annotations & Data Sources

| File | Source |
|------|--------|
| `banc_meta.feather` | `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/` |
| FAFB annotations | [flyconnectome/flywire_annotations](https://github.com/flyconnectome/flywire_annotations) |
| MANC annotations | [male-cns.janelia.org/download](https://male-cns.janelia.org/download/) |

---

## References

- Schlegel et al. (2024) *Nature* [doi:10.1038/s41586-024-07686-5](https://doi.org/10.1038/s41586-024-07686-5)
- Dorkenwald et al. (2024) *Nature* [doi:10.1038/s41586-024-07558-y](https://doi.org/10.1038/s41586-024-07558-y)
- Pospisil et al. (2024) *Nature* [doi:10.1038/s41586-024-07982-0](https://doi.org/10.1038/s41586-024-07982-0)
- Bates et al. (2025) *bioRxiv* [doi:10.1101/2025.07.31.667571](https://doi.org/10.1101/2025.07.31.667571)
- Takemura et al. (2024) *eLife* [doi:10.7554/eLife.97769](https://doi.org/10.7554/eLife.97769)
- Witvliet et al. (2021) *Nature* [doi:10.1038/s41586-021-03778-8](https://doi.org/10.1038/s41586-021-03778-8)
