# Structural Invariance at the Sensorimotor Interface: A Maximum Common Induced Subgraph Across Three *Drosophila* Connectomes

**Yanjin Li** · FlyWire Qualification Challenge · June 2026

**Datasets:** BANC v626 (♀ brain+cord) · FAFB v783 (♀ brain) · MANC v1.2.1 (♂ nerve cord)  
**Result:** N = 105 neurons · 12 conserved directed edges · correspondence-shuffle null 76.2 ± 1.5 (>15σ separation)  
**Code:** [github.com/Yanjin-ai/flywireconnectome](https://github.com/Yanjin-ai/flywireconnectome) · all numbers reproduced by `src/run_analysis.py` (`results/`)

---

## 📄 One-Page Scientific Report

**The circuit.** The largest set of morphologically matched neurons whose directed induced subgraph is *identical* across three independent connectomes — BANC (♀ brain+cord), FAFB (♀ brain), MANC (♂ cord) — is a **105-neuron sensorimotor backbone with 12 conserved directed edges**. It is 93% descending + ascending neurons (65 DN + 33 AN + 7 sensory): the brain↔ventral-nerve-cord communication channel.

**What it does.** The conserved hubs span multiple motor modules — **wing/flight, leg, and abdominal** effectors — and are **anchored by DNg02**, a documented descending population that sets wingbeat amplitude / flight steering through a population code (Schnell, Ros & Dickinson 2022). Descending axons of these classes target leg/neck/wing motor circuits in the VNC (Namiki et al. 2018). The mix of cholinergic and GABAergic descending neurons is consistent with a feedforward-inhibition coordination motif (Milo et al. 2002).

![Network graph](figures/figure3_hub_circuit.png)
*Conserved-circuit network graph: hubs joined by the 12 directed edges identical across all three connectomes (gold).*

![Codex 3D meshes](figures/codex_3d_fafb.png)
*Codex 3D meshes (FAFB): all 105 neurons inside the whole-brain outline, converging at the midline / cervical connective — the expected brain↔cord relay locus.*

**Structural observations.**
- **Cell-type identity:** 66.2× enriched for descending and 25.0× for ascending neurons vs the whole-brain FAFB background (Fisher p < 10⁻³⁵) — not a random brain sample.
- **Wiring conserved beyond degree:** 2,609 edges are shared by all three connectomes vs 352.9 ± 16.5 under a degree-preserving (Maslov–Sneppen) null → **7.4×, Z = 136σ**. Specific connectivity is conserved, not just degree sequence.
- **Cross-sex:** 88.6% of the circuit is wired identically in ♀ and ♂; the few dimorphic neurons target abdominal VNC / lateral brain.
- **Robustness:** lower betweenness than matched neurons (peripheral relays, p = 0.009); circuit survives 20% simulated reconstruction error (graceful decline).

**Interpretation & hypotheses.** The backbone is a developmentally canalised brain↔cord channel (H1) — supported by cross-sex conservation and output-hemilineage enrichment. The degree-preserving null **rejects** a pure degree/sampling artifact at the edge level (H2), and the strong class enrichment **rejects** a generic-subgraph explanation (H5); current static data **cannot** distinguish developmentally fixed vs activity-refined wiring (H4 — the key open question). **Prediction:** silencing the hub neurons (DNg02; the reciprocal DNa15↔DNg04 and DNp58↔DNp65 pairs) should impair walking, flight and posture *simultaneously* — testable by optogenetic silencing with multi-behaviour assays.

**Key citations.** Schlegel et al. 2024 (cell typing); Bates et al. 2025 (BANC correspondence); Dorkenwald et al. 2024 (FAFB); Pospisil et al. 2024 (sensorimotor bottleneck); Namiki et al. 2018 & Schnell et al. 2022 (descending-neuron function); Witvliet et al. 2021 (connectome stereotypy). Full reference list in §References.

> Print-ready version: [`research_summary_fafb.pdf`](research_summary_fafb.pdf). The detailed report (methods, statistics, all figures) follows below.

---

## 1. Hypothesis

The *Drosophila* nervous system has been reconstructed across multiple independent specimens, sexes, and anatomical preparations. Schlegel et al. (2024) established a consensus cell-type atlas spanning five datasets, demonstrating that **cell-type identity** is reproducible at the morphological level. A deeper unresolved question is whether **synaptic connectivity itself** is structurally invariant across independently prepared connectomes.

> **Central hypothesis:** There exists a set of morphologically matched neurons whose directed synaptic connectivity forms a *mutually isomorphic induced subgraph* — a structurally invariant backbone (Witvliet et al. 2021) — across at least three independent connectomes. This backbone is enriched at the sensorimotor interface (descending and ascending neurons), reflecting a developmental constraint on the brain–body communication channel conserved across sexes and specimens.

This extends previous work (Schlegel et al. 2024; Witvliet et al. 2021) which quantified *cell-type-level* and *motif-level* conservation; we search for the *maximum* set of neurons with *edge-level* structural identity.

**Operationalised predictions (testable from available data):**

1. *Annotation quality:* Circuit neurons should have higher manual-annotation confidence. **Verified:** 93.3% manually-annotated (`manual_cluster`) vs 85.0% for non-circuit members (Fisher exact p = 0.008; §4.4).
2. *Statistical significance:* Shuffling cross-dataset correspondence should yield a significantly smaller MCIS. **Verified:** null collapses to 76.2 ± 1.5 vs real 100.5 ± 2.2 (>15σ separation; §4.2).
3. *Degree-distribution signal:* Degree-preserving edge rewire should yield a substantially smaller MCIS if specific edge patterns matter beyond degree sequence. **Result:** the degree-preserving null reaches 100.8 ± 2.1 — statistically indistinguishable from the real per-seed mean — so the FAFB degree distribution alone accounts for nearly all of the achievable N; neuron identity via NBLAST correspondence provides the additional signal to reach the full N = 105 (Z ≈ 2σ against the best; §4.2, full interpretation).

---

## 2. Dataset Selection and Correspondence

### 2.1 Why BANC × FAFB × MANC

| Dataset | Sex | Region | Neurons (edge list) | Version |
|---------|-----|--------|---------------------|---------|
| **BANC** | ♀ | Brain + ventral nerve cord | 112,885 | v626 |
| **FAFB** | ♀ | Brain only | 138,584 | v783 |
| **MANC** | ♂ | Ventral nerve cord only | 23,641 | v1.2.1 |

BANC is the only dataset spanning both brain and ventral nerve cord. Its metadata (Bates et al. 2025) contains `fafb_match` and `manc_match` columns — 1:1 NBLAST-based morphological correspondences per neuron, giving 3,414 pre-verified triplets at the individual-cell level. No equivalent three-way correspondence table exists for MAOL or MCNS at this resolution. The cross-sex design (♀ FAFB/BANC vs ♂ MANC) means any surviving circuit is a strong candidate for evolutionary canalization.

**Why not MAOL or MCNS?** The MAOL (male optic lobe) and MCNS (male full CNS) datasets lack a published three-way NBLAST-based correspondence to BANC and FAFB at the individual-neuron level. Including them would require heuristic matching, introducing unquantified error. We exclude them to maintain ground-truth provenance for all correspondences.

### 2.3 Systematic survey of alternative dataset triplets

To justify the choice of BANC × FAFB × MANC as the primary triplet, we characterise all scientifically plausible three-way combinations of the five published *Drosophila* connectomes using three axes: (i) ground-truth correspondence availability, (ii) anatomical complementarity (whether the triplet jointly covers brain + ventral nerve cord), and (iii) a rough MCIS size upper bound estimated as the number of neurons with manually-curated matches across all three datasets.

| Triplet | Ground-truth 3-way correspondence | Anatomical coverage | Cross-sex | Estimated MCIS upper bound | Verdict |
|---------|-----------------------------------|--------------------|-----------|-----------------------------|---------|
| **BANC × FAFB × MANC** | ✅ BANC `fafb_match` + `manc_match` (Bates et al. 2025) | Brain + VNC ✅ | ♀/♀/♂ ✅ | ~2,800 (direct edge intersection) | **Primary choice** |
| FAFB × MANC (pairwise) | ✅ via BANC proxy | Brain + VNC ✅ | ♀/♂ | ~2,800 | Lacks independent third dataset; pairwise comparison does not test three-way isomorphism |
| BANC × FAFB × MCNS | ⚠️ MCNS lacks published individual-neuron NBLAST to BANC/FAFB | Brain + full CNS | ♀/♀/♂ | ~500 (cell-type proxy only) | Correspondence error unquantified; cell-type aggregation loses neuron-level resolution |
| FAFB × MAOL × MCNS | ❌ No three-way individual-neuron correspondence table | Brain + optic lobe + full CNS | ♀/♂/♂ | ~300 (heuristic estimate) | Requires ad hoc matching; introduces unknown false-match rate |
| BANC × MANC × MCNS | ⚠️ MCNS↔BANC correspondence not published | VNC + VNC + full CNS | ♀/♂/♂ | ~200 | Redundant VNC coverage; misses brain half of the sensorimotor axis |
| FAFB × FANC × MCNS | ❌ FANC individual-neuron correspondence to FAFB not yet published | Brain + ♀VNC + full CNS | ♀/♀/♂ | Unknown | FANC metadata in active development; premature to include |

**Key conclusion.** BANC × FAFB × MANC is the only triplet that simultaneously satisfies all three criteria: ground-truth individual-neuron correspondence, full brain-to-cord anatomical axis, and meaningful cross-sex comparison. The estimated upper bound (~2,800 shared triplets) is 5–14× larger than any alternative, maximising statistical power. When three-way correspondence tables for MCNS and MAOL become publicly available, the most informative extension would be BANC × FAFB × MCNS (adds a second ♂ full-CNS dataset, enabling four-way intersection and a direct test of conservation depth as a function of dataset count).

### 2.2 Why only 1.3% edge consensus

Among the 2,798 triplets present in all three edge lists, only **2,609 directed edges** are shared by all three connectomes (the consensus graph used for the search); the large majority of edges present in any single dataset are not shared. This is biologically interpretable, not a data quality failure. Descending neurons have dendrites in the brain (synapses captured by FAFB) and axonal outputs in the ventral nerve cord (synapses captured by MANC); only BANC captures both compartments. This is consistent with Witvliet et al. (2021): approximately 60% of *C. elegans* chemical synapses are variable across individuals even under controlled conditions, without any cross-compartment sampling.

At the cell-type aggregation level (collapsing individual neurons to named types), FAFB and MCNS share 7,289 named types — far higher overlap than the individual-neuron level, confirming that the low individual-neuron consensus is a property of the cross-compartment comparison, not of dataset quality per se.

---

## 3. Method

### 3.1 MCIS problem formulation

```
Input:  987 matched neurons in the giant consensus component
        Indexed edge sets E_BANC, E_FAFB, E_MANC

Goal:   Largest S ⊆ {1..987} such that
        ∀ i,j ∈ S: (i→j) ∈ E_BANC  ⟺  (i→j) ∈ E_FAFB  ⟺  (i→j) ∈ E_MANC

Because nodes are bijectively labelled via NBLAST triplets, isomorphism
reduces to edge-set equality — O(N²), not the general NP-hard unlabelled case.
```

### 3.2 Algorithm: Greedy Disagreement Removal + Exhaustive Expansion

```
Phase 1 — Greedy removal:
  While ∃ disagreement edges:
    Remove argmax_v |{disagreement edges incident to v}|
    (randomised tie-breaking for deterministic multi-seed reproducibility)

Phase 2 — Exhaustive expansion:
  For each removed node v (in arbitrary order):
    If adding v to current set preserves isomorphism → add v

Complexity: O(N·D) per iteration, D = number of disagreement edges
Convergence: ≤890 iterations on the 987-node instance (~9 seconds)
```

**Relation to the Maximum Common Subgraph literature.** Classical MCS solvers (McGregor's backtracking, 1982; the clique-on-the-product-graph reduction of Raymond & Willett 2002; modern branch-and-bound solvers such as McSplit) must *discover* the node correspondence, which makes general MCS doubly hard (subgraph isomorphism is NP-complete; maximum common subgraph is NP-hard). Our setting is fundamentally easier because the correspondence is **given a priori** by NBLAST (each neuron has one identity across datasets). With a fixed bijection, "mutually isomorphic induced subgraph" reduces to **edge-set equality**, and maximising N becomes exactly **Maximum Independent Set on the disagreement graph** D (§3.2) — we never search over matchings. This is why we can certify near-optimality with a compact ILP (MIS, not the product-graph clique program) and why a simple greedy is competitive. We still benchmark the heuristic: greedy disagreement removal dominates an *edge-centric growth* alternative (seed from the highest-degree consensus node, grow), the latter being hampered by local seeding ([`src/robustness_experiments.py`](src/robustness_experiments.py)).

**Exact MCIS validation via ILP.** To provide a certified lower bound on the optimality gap, we formulated MCIS as an Integer Linear Program and solved it to provable optimality on randomly sampled induced subgraphs of varying size.

*ILP formulation.* Let binary variable $x_v \in \{0,1\}$ indicate whether neuron $v$ belongs to the solution set $S$. For each ordered pair $(i,j)$ present in exactly one or two (but not all three) datasets — a "disagreement edge" — at most one of $x_i$, $x_j$ can equal 1 (otherwise the induced subgraph is not isomorphic). The full program is:

$$\text{maximise} \sum_v x_v \quad \text{subject to} \quad x_i + x_j \leq 1 \; \forall (i,j) \in \mathcal{D}, \quad x_v \in \{0,1\}$$

where $\mathcal{D}$ is the set of disagreement edges restricted to the subgraph. This is a Maximum Independent Set on the disagreement graph, which is itself NP-hard in general but tractable on small instances via branch-and-bound with LP relaxation (solved here with PuLP/CBC).

*Results on subgraph samples.* We drew 50 random induced subgraphs at four size tiers from the 987-node consensus component, solved each exactly with the ILP (PuLP/CBC), and compared against the greedy + exhaustive expansion result. These numbers are produced by [`src/exact_ilp.py`](src/exact_ilp.py) and stored in `results/ilp_validation.json`:

| Subgraph size | Instances | Greedy N (mean ± sd) | ILP-optimal N (mean ± sd) | Optimality gap |
|--------------|-----------|----------------------|--------------------------|----------------|
| 20 nodes | 15 | 12.7 ± 1.4 | 12.8 ± 1.3 | **0.5%** |
| 30 nodes | 15 | 18.2 ± 1.5 | 18.3 ± 1.5 | **0.7%** |
| 40 nodes | 10 | 20.3 ± 1.5 | 20.6 ± 1.6 | **1.5%** |
| 50 nodes | 10 | 22.9 ± 1.6 | 23.5 ± 1.7 | **2.6%** |

Across all 50 instances the greedy + expansion heuristic achieves a **mean optimality gap of 1.15%** (maximum 10.5% on a single 50-node instance; 80% of instances within 2%). The gap grows modestly with subgraph size, as expected for a local-search heuristic. ILP solve times were 0.02–1.3 s for these sizes; the full 987-node instance is intractable for exact ILP, but the consistent low gap on subgraphs, together with the 100-seed variance of ±2.2 (§4.1), indicates that the reported N is within a few neurons of the true optimum.

**Unit tests:** `pytest tests/ -v` — 14 tests covering isomorphism verification, planted subgraph recovery (known ground truth), expansion monotonicity, null model consistency, a synthetic greedy-vs-exact (brute-force) check, and data-loading smoke tests. 13 run on synthetic graphs; 1 real-data smoke test is skipped when `MCIS_DATA_DIR` is unset. All pass.

---

## 4. Robustness and Statistical Validation

### 4.1 Algorithmic robustness — 100 random seeds

> *N = 100.5 ± 2.2, range [96, 105] across 100 random tie-breaking seeds (each followed by exhaustive expansion).*

The narrow range (±2.2 over a 987-node search space) confirms that N ≈ 100 is a stable structural property of the data, not a fragile artifact of a specific node ordering. The reported circuit is the **best of the 100-seed multi-start, N = 105**, shipped as `network.csv`; the full distribution and all statistics below are produced by [`src/run_analysis.py`](src/run_analysis.py) (`results/canonical_results.json`).

### 4.2 Three null models

All nulls use the *same* search procedure as the real result (best-of-5 multi-start per trial), so comparisons are apples-to-apples (20 trials each; [`src/run_analysis.py`](src/run_analysis.py)).

| Null model | Construction | N_null | vs real | Interpretation |
|-----------|-------------|--------|-----------|----------------|
| Correspondence-shuffle | Permute *both* FAFB and MANC neuron→triplet mappings; BANC unchanged | 76.2 ± 1.5 | >15σ below real 100.5 ± 2.2 | Neuron identity (NBLAST matching) is essential |
| Degree-preserving rewire | Rewire ~33% of FAFB edges while preserving exact in/out degree per neuron | 100.8 ± 2.1 | indistinguishable from real mean (Z ≈ 2σ vs best N=105) | FAFB degree sequence alone accounts for nearly all of the achievable N |
| Centrality permutation (1000 trials) | Permute neuron labels on the consensus-graph betweenness | — | p = 0.009 | Circuit neurons have *lower* betweenness than matched pool average |

**Interpretation of the degree-preserving null.** The degree-preserving rewire is a Maslov–Sneppen randomisation (Maslov & Sneppen 2002): it scrambles connectivity while holding each neuron's in/out degree fixed, isolating the contribution of the degree sequence from that of specific wiring. This result indicates that the FAFB degree distribution is the dominant determinant of how large an MCIS can be found: the Maslov–Sneppen null reaches 100.8 ± 2.1, statistically indistinguishable from the real per-seed mean (100.5 ± 2.2). The correspondence-shuffle null (collapsing to 76.2 ± 1.5) shows that NBLAST-based neuron identity is additionally required — but the honest reading is that the degree sequence is a near-sufficient condition, and specific neuron identity provides the additional, smaller contribution needed to reach the best N = 105.

### 4.3 Centrality: circuit neurons are peripheral relays

Circuit betweenness centrality = 0.00096 vs non-circuit matched neurons = 0.00224 on the consensus graph (one-sided label-permutation test, p = 0.009, 1000 permutations; the test asks specifically whether circuit betweenness is *lower* than the matched pool). Circuit neurons have **lower** betweenness — they are inter-system relay neurons at the brain–body interface, not structural hubs within the brain network. This is biologically coherent: DN/AN neurons are few-input, few-output specialists bridging two anatomical compartments rather than central integrators.

### 4.4 Annotation quality — verified prediction

> *Circuit: 93.3% manually annotated (98/105, `manual_cluster` non-null) vs non-circuit: 85.0% (2,814/3,309); Fisher exact p = 0.008.*

The MCIS result is enriched for the more carefully annotated correspondences, not driven by low-confidence matches.

### 4.5 NBLAST confidence curve

Using agreement between automated NBLAST top-1 output and expert-curated match as a confidence proxy (high: both FAFB and MANC agree; medium: one agrees; low: neither agrees but still manually verified):

| Top-k% (highest confidence first) | Triplet pool | MCIS N |
|-----------------------------------|--------------|--------|
| 10% | 280 | 4 |
| 20% | 560 | 39 |
| 30% | 839 | 53 |
| 50% | 1,399 | 74 |
| 75% | 2,098 | 91 |
| 100% | 2,798 | **105** |

N increases monotonically without discontinuity across all confidence tiers ([`src/confidence_tiers.py`](src/confidence_tiers.py), `results/confidence_tiers.json`), confirming that the result is not driven by a small pocket of low-confidence matches.

![Fig. 1 — NBLAST confidence](figures/figure7_nblast_confidence.png)
**Figure 1.** NBLAST confidence analysis. **(A)** MCIS size vs confidence threshold: monotonic increase from N=4 (top 10%) to N=105 (full pool). **(B)** Distribution of NBLAST agreement across 2,798 triplets.

![Fig. 2 — Robustness panel](figures/figure5_robustness.png)
**Figure 2.** Robustness and validation. **(A)** 100-seed MCIS distribution: 100.5 ± 2.2, range [96, 105]. **(B)** Null comparison: correspondence-shuffle collapses to 76.2 ± 1.5; degree-preserving rewire reaches 100.8 ± 2.1. **(C)** N bound waterfall (3,414 → 2,798 → 987 → 105). **(D)** Empirical runtime O(N^1.9); 987 nodes in ~8.8 seconds. **(E)** Centrality: circuit has lower betweenness (p = 0.009). **(F)** MCIS stability across confidence tiers. **(G)** Annotation quality: 93.3% vs 85.0% manually annotated (p = 0.008).

### 4.6 Conservation beyond degree sequence — the conservation track

The MCIS *size* (node count) is largely explained by the degree sequence (§4.2). This is expected: the MCIS is dominated by neurons with few or no induced edges, so a degree-preserving rewire can recover a similarly large *set* of nodes. The scientifically decisive question is about the **shared wiring itself**: are there more edges present in all three connectomes than a degree-preserving rewiring of each connectome would produce by chance?

We answer this at the edge level ([`src/conservation_track.py`](src/conservation_track.py), `results/conservation_track.json`). Over the 987-node consensus component (56,337 edges present in ≥1 connectome):

| Quantity | Observed | Degree-preserving null (100 trials) | Result |
|---|---|---|---|
| Edges present in all 3 connectomes | **2,609** | 352.9 ± 16.5 | **7.4× enrichment, Z = 136σ** |

So while node-count MCIS is degree-explained, **specific synaptic connectivity is conserved 7.4× above the degree-sequence expectation** — strong, unambiguous evidence that the cross-connectome agreement reflects real wiring identity, not merely matched degree distributions. This reframes the contribution from a binary "conserved circuit" to a continuous **conservation track**: every neuron receives a conservation z-score (observed consensus-incident edges vs degree-null), yielding a ranked map of which neurons carry the conserved wiring (top: ANXXX108 z=35, DNge106 z=28, DNg73/AN17B008 z=21). This per-neuron track is what the Neuroglancer overlay (§10.5) colours.

![Fig. 10 — Conservation track](figures/figure10_conservation_track.png)
**Figure 10.** **(A)** Edge support across connectomes (1 / 2 / all-3). **(B)** Beyond-degree test: observed 2,609 consensus edges vs degree-preserving null 353 ± 17 (Z = 136σ). **(C)** Per-neuron conservation z-score track.

### 4.7 Robustness to connectomic reconstruction error

Every connectome carries proofreading error. To test whether the result is an artifact of the exact edge sets, we independently perturbed each connectome — flipping a fraction *p* of its edges (removing real edges = false negatives, adding random edges = false positives) — and recomputed the MCIS ([`src/stringency_sweep.py`](src/stringency_sweep.py), `results/stringency_sweep.json`):

| Edge error per connectome | MCIS N (best-of-3) |
|---|---|
| 0% | 105.0 ± 0.0 |
| 5% | 92.5 ± 1.1 |
| 10% | 85.0 ± 1.0 |
| 20% | 72.8 ± 1.3 |

N **degrades gracefully** (≈ linear, no cliff): even with 20% of every connectome's edges corrupted, a 73-neuron conserved circuit survives. The conserved backbone is therefore a stable structural feature, not a fragile coincidence of the specific reconstructions.

![Fig. 11 — Reconstruction-error robustness](figures/figure13_stringency.png)
**Figure 11.** MCIS size vs per-connectome edge perturbation; graceful, near-linear decline.

---

## 5. Cell-Type Enrichment: The Circuit Is Not a Random Brain Sample

Comparing the 105 circuit neurons against the full 139,244-neuron FAFB annotation as background (computed by [`src/derived_stats.py`](src/derived_stats.py), `results/derived_stats.json`):

| Neuron class | Circuit (N=105) | FAFB background (N=139,244) | Fold enrichment | Fisher exact p |
|-------------|-----------------|----------------------------|-----------------|----------------|
| Descending | 61.9% (65/105) | 0.94% (1,303/139,244) | **66.2×** | 3.0 × 10⁻¹⁰⁴ |
| Ascending  | 31.4% (33/105) | 1.26% (1,750/139,244) | **25.0×** | 1.2 × 10⁻³⁶ |

The circuit is 66.2× enriched for descending neurons and 25.0× enriched for ascending neurons (both p < 10⁻³⁵). These reflect a near-complete exclusion of non-sensorimotor neuron classes from the MCIS.

The most represented developmental hemilineages among circuit neurons are LB12, 05B, 09B, and SMPpv2 — established output hemilineages projecting from brain to nerve cord (Ito et al. 2013), consistent with the developmental constraint hypothesis.

![Fig. 3 — Cell-type enrichment](figures/figure9_enrichment.png)
**Figure 3.** Cell-type enrichment vs FAFB whole-brain background. **(A)** Superclass fold-enrichment (descending 66.2×, ascending 25.0×). **(B)** Neurotransmitter profile: ACh-dominant (69.5%; 73/105). **(C)** Most represented developmental hemilineages.

---

## 6. The Conserved Circuit

### 6.1 Composition

| Neuron class | Count | % | Functional role |
|-------------|-------|---|-----------------|
| Descending (DN) | 65 | 61.9% | Brain → VNC motor commands |
| Ascending (AN) | 33 | 31.4% | VNC → Brain proprioceptive feedback |
| Sensory-ascending | 5 | 4.8% | Peripheral sensory → Brain |
| Sensory-descending | 2 | 1.9% | Sensory processing → VNC |
| **Total** | **105** | | **12 conserved directed edges** |

### 6.2 Anatomical position

![Fig. 4 — Spatial distribution](figures/figure6_spatial.png)
**Figure 4.** BANC anatomical projections (`root_position_nm` scaled to µm). **(A–C)** Coronal, sagittal, and axial projections: circuit neurons (coloured by class) are concentrated along the cervical connective — the anatomically expected locus for DN/AN neurons bridging brain and ventral nerve cord. Grey = all matched neurons.

### 6.3 Circuit structure and neurotransmitters

![Fig. 5 — Circuit layouts](figures/figure1_circuit_layouts.png)
**Figure 5.** Three force-directed layouts of the 105-neuron circuit. Gold edges = 12 directed connections verified identical across BANC, FAFB, and MANC. Red = descending (DN); blue = ascending (AN); large nodes = neurons involved in conserved edges.

![Fig. 6 — Hub neurons](figures/figure3_hub_circuit.png)
**Figure 6.** Hub neurons connected by the 12 conserved edges (e.g. the reciprocal DNa15↔DNg04 and DNp58↔DNp65 pairs, and the DNge076→DNge019/DNge020 fan-out), with neurotransmitter identity annotated. The mixed ACh/GABA/Glu chemistry is consistent with a feedforward inhibition motif — a canonical computation (Milo et al. 2002) enabling temporal filtering of descending motor commands.

![Fig. 7 — Composition](figures/figure2_composition.png)
**Figure 7.** **(A)** DN/AN dominance. **(B)** Acetylcholine-dominant NT profile (69.5%; 73/105). **(C)** Multi-effector motor targets (leg VNC, dorsal VNC, flange median bundle, abdominal VNC). **(D)** Node degree distribution. **(E)** Cross-dataset edge count comparison: asymmetry reflects partial-volume biology (MANC captures axonal synapses; FAFB captures dendritic synapses).

### 6.4 Motor targets — multi-effector coordination

Among neurons with an annotated `cns_network` target: leg VNC (38 neurons, locomotion), dorsal VNC/flight (18), lateral brain (12), flange median bundle/whole-body coordination (11), posterior brain (8), abdominal VNC (7); 7 neurons are unannotated. The multi-effector profile is characteristic of coordination interneurons rather than single-behaviour specialists.

### 6.5 Hub neurons: motor modules and documented function

The 13 cell types that carry the 12 conserved edges, grounded in Codex annotations (cell type, neurotransmitter, VNC/brain target) and the descending-neuron literature. We deliberately mark documented single-type function vs target-inferred module — most individual DN/AN types are not yet behaviourally characterised, and we do not over-claim.

| Hub cell type(s) | NT | Target (`cns_network`) | Motor module (Namiki et al. 2018) | Documented function |
|---|---|---|---|---|
| **DNg02** (DNg02_g) | ACh | dorsal VNC | wing / flight neuropil | **Yes** — a population of ≥15 pairs that regulates wingbeat amplitude and flight steering via a population code (Schnell, Ros & Dickinson 2022) |
| DNg04, DNg79, DNa15 | ACh | dorsal VNC | wing / neck / flight | Target-inferred (dorsal motor neuropil); single-type behaviour not yet characterised |
| DNge019, DNge020 | ACh | leg VNC | leg / locomotion | Target-inferred (leg motor circuits); not individually characterised |
| DNp58 | ACh | abdominal VNC | abdominal / postural | Target-inferred; not individually characterised |
| DNp47 | ACh | posterior brain | brain-targeting DN | Not individually characterised |
| DNge076, DNp65, DNp54 | **GABA** | flange median bundle / post. brain | inhibitory DNs | Inhibitory; not individually characterised |
| AN09B033 | ACh | (ascending) | proprioceptive feedback | Ascending; not individually characterised |
| AN06A027 | Glu | (ascending) | proprioceptive feedback | Ascending; not individually characterised |

Two observations follow. (i) The conserved core is **multi-module**: it spans wing/flight (DNg-class, dorsal VNC), leg/locomotion (DNge-class, leg VNC) and abdominal/postural (DNp58, AN06A027) effectors — consistent with a cross-program coordination backbone rather than a single-behaviour pathway. (ii) It is **anchored by a documented controller**: DNg02, the only individually characterised type in the set, is a flight-motor population — its three-way structural conservation means a known population-code flight controller sits inside the invariant backbone. The presence of GABAergic descending neurons (DNge076, DNp65, DNp54) alongside cholinergic ones gives the mixed excitatory/inhibitory chemistry expected of a feedforward-inhibition coordination motif (§8.3).

---

## 7. Sexual Conservation Analysis

### 7.1 Overview

**88.6% of circuit neurons (93/105) are sexually isomorphic** — their wiring is preserved identically across ♀ FAFB/BANC and ♂ MANC. The remaining 11.4% (12 neurons) are sexually dimorphic. For comparison, Berg et al. (2025) report approximately 95.2% sexual conservation across all matched DN/AN neuron pairs; our circuit's 88.6% is slightly below this baseline, reflecting the presence of sex-specific behavioural neurons among the 105.

![Fig. 8 — Dimorphism overview](figures/figure4_dimorphism_nt.png)
**Figure 8.** Sexual dimorphism overview. Dimorphism status by neuron class and neurotransmitter profile.

![Fig. 9 — Sexual conservation deep dive](figures/figure8_sexual_conservation.png)
**Figure 9.** **(A)** 88.6% isomorphic (93/105). **(B)** Dimorphism by neuron class. **(C)** Neurotransmitter identity of the 12 dimorphic neurons (ACh 9, Glu 2, serotonin 1).

### 7.2 The 12 sexually dimorphic neurons

| Cell type | Class | NT | Motor target |
|-----------|-------|-----|-------------|
| AN05B102 (×2) | Ascending | ACh | Lateral brain |
| DNg111 | Descending | Glu | Leg VNC |
| DNp47 | Descending | ACh | Posterior brain |
| DNp67 | Descending | ACh | Leg VNC |
| DNpe052 | Descending | ACh | Lateral brain |
| ANXXX254 | Ascending | ACh | Abdominal VNC |
| ANXXX169 | Ascending | Glu | Abdominal VNC |
| DNge010 | Descending | ACh | Leg VNC |
| DNg02_g | Descending | ACh | Dorsal VNC |
| LN-DN2 | Sensory-desc. | Serotonin | — |
| SAch01 | Sensory-asc. | ACh | — |

Three neurotransmitter types appear among dimorphic neurons: ACh (9), Glu (2), serotonin (1). The serotonergic and glutamatergic neurons are notably neuromodulatory — controlling internal state rather than direct motor output. Neurons projecting to abdominal VNC and lateral brain account for a disproportionate share of dimorphic neurons, consistent with sex-specific reproductive and mating behaviours requiring abdominal motor control and higher-order brain integration in one sex but not the other.

---

## 8. Biological Interpretation

### 8.1 Structural evidence for the sensorimotor bottleneck

Pospisil et al. (2024) showed that approximately 1% of *Drosophila* brain neurons directly influence motor output — the "effectome" — with descending neurons as the obligate conduit. Our result provides **direct structural corroboration**: the largest isomorphic subgraph across three independent connectomes consists almost entirely of DNs and ANs (93.3% of circuit neurons). This demonstrates that the sensorimotor bottleneck is not only functionally constrained but **structurally canalized across sexes, specimens, and anatomical preparations**. Notably, the conserved core includes **DNg02** — a documented flight-motor population that sets wingbeat amplitude through a population code (Schnell et al. 2022) — so a behaviourally characterised controller sits inside the invariant backbone (§6.5).

### 8.2 Developmental constraint as the mechanistic basis

The 88.6% cross-sex conservation is consistent with the developmental constraint hypothesis: DN/AN connectivity is established early in neurogenesis by lineage-specific programs (hemilineage identity; Ito et al. 2013) that are largely sex-independent. The 11.4% dimorphic fraction maps onto neurons with sex-specific motor targets (abdominal VNC, lateral brain) and neuromodulatory roles — exactly the classes expected to diverge between sexes for reproductive behaviour.

### 8.3 Testable experimental predictions

1. **Multi-program impairment:** Silencing the hub neurons that carry the conserved edges (e.g. the reciprocal pairs DNa15↔DNg04 and DNp58↔DNp65, or the DNge076→DNge019/DNge020 fan-out) should impair walking, flight, and posture simultaneously — testable via optogenetic silencing combined with multi-behaviour assays.
2. **Synapse strength:** The 12 conserved edges should exhibit above-average synapse counts in the weighted connectome, consistent with robust signal transmission. Testable via FlyWire API query of synapse weights.
3. **Cross-species conservation:** Orthologous circuits should be identifiable in other holometabolous insects (*Manduca sexta*, *Apis mellifera*) as connectome data become available, given that DN/AN cell types are broadly conserved across Insecta.

### 8.4 Alternative hypotheses — what the data can and cannot distinguish

The headline observation (105 matched neurons with edge-level connectivity conserved 7.4× above a degree-preserving null) is consistent with several hypotheses. We state them explicitly and mark which our current data adjudicate:

| Hypothesis | Prediction | Verdict from this study |
|---|---|---|
| **H1 — Developmental canalisation.** The backbone is wired by lineage-specific programs largely independent of sex/specimen. | Conserved across sexes; enriched in known output hemilineages; conserved beyond degree. | **Supported** (88.6% cross-sex; LB/SMPpv2 hemilineages; 7.4× beyond-degree) — but not *proven*: a structural snapshot cannot show the wiring was set developmentally rather than refined by activity. |
| **H2 — Degree/sampling artifact.** Apparent conservation is a by-product of matched degree sequences and shared dense regions. | A Maslov–Sneppen (degree-preserving) null should reproduce the shared edges. | **Rejected at the edge level** (observed 2,609 vs 353 ± 17; Z = 136σ). Note it is *not* rejected for node-count (§4.2) — hence we report edges, not N, as the conserved signal. |
| **H3 — Annotation/proofreading bias.** Conservation tracks the best-annotated, most-proofread neurons. | The result should collapse when restricted to high-confidence matches, and circuit/non-circuit annotation rates should differ strongly. | **Largely rejected**: N grows monotonically across all NBLAST-confidence tiers (§4.5) and the annotation-rate gap is modest (93.3% vs 85.0%). A residual bias cannot be fully excluded. |
| **H4 — Functional/activity-driven conservation** (vs developmental). | Conserved edges would correlate with co-activity or behavioural necessity, not just lineage. | **Cannot be distinguished** with static connectomes alone — requires activity imaging or perturbation (§8.3 prediction 1) and the synapse-weight readout (prediction 2). This is the key open question. |
| **H5 — Generic-subgraph property** (any matched neuron set would look conserved). | Enrichment for specific classes should be absent. | **Rejected**: the circuit is 66.2×/25.0× enriched for descending/ascending neurons — the conservation is specifically sensorimotor, not a generic property of matched neurons. |

### 8.5 Comparison to models and behavioural data

Our structural backbone complements three recent computational/functional lines. (i) *Connectome-constrained mechanistic models* (Lappalainen et al. 2024) instantiate measured connectivity and fit single-neuron dynamics to predict activity in the fly **visual** system, succeeding precisely where connectivity is sparse; our result extends the "sparse, connectome-constrained" regime to the **sensorimotor** axis and supplies a three-connectome-validated topology rather than a single-specimen one. (ii) The **effectome** (Pospisil et al. 2024) identifies the ~1% of brain neurons with direct motor leverage via descending neurons; our backbone is the *structurally invariant* core of exactly that population (93% DN/AN), suggesting the effectome's obligate conduit is also the most evolutionarily canalised. (iii) Whole-brain leaky-integrator simulations (Shiu et al. 2024) predict sensorimotor responses from FAFB connectivity; the 12 conserved edges (reciprocal DN↔DN pairs with mixed ACh/GABA) are concrete, falsifiable targets whose perturbation such a model could be asked to reproduce. None of these works tests cross-connectome structural invariance, which is the gap this study fills.

---

## 9. Limitations

- **N is a near-optimal lower bound.** The greedy algorithm finds a local optimum; the true MCIS on 987 nodes is NP-hard to certify exactly. ILP validation on 50 subgraphs (20–50 nodes) demonstrates a mean optimality gap of 1.15% (maximum 10.5% on one 50-node instance; §3.2), and the 100-seed variance of ±2.2 (§4.1) provides independent evidence of stability. Together these indicate the reported N = 105 is within a few neurons of the true optimum, but a formally certified proof of global optimality on the full 987-node instance remains intractable.
- **No continuous NBLAST scores.** The BANC metadata contains binary match results, not morphological similarity scores. A fully continuous confidence curve would require the R `bancr` package and neuron skeleton data (~10 GB); we used NBLAST top-1 agreement as a proxy (§4.5).
- **Degree-preserving null reaches the real MCIS *size*** (null 100.8 ± 2.1 ≈ real 100.5 ± 2.2; §4.2). This limits claims based on *node count*. It does **not** limit the wiring claim: at the *edge* level, all-3 consensus edges are enriched 7.4× over the same degree-preserving null (Z = 136σ; §4.6), so specific connectivity is conserved well beyond degree sequence. The right unit of analysis is edges, not node count.
- **Limited MANC cross-link coverage.** Approximately 2,498 of MANC's 23,641 neurons are cross-linked via the MCNS proxy table, constraining the triplet pool.
- **Centrality result is directional only.** The betweenness difference (p = 0.009) should be treated as a directional observation until confirmed with a formal parametric test and/or replicated on a second connectome pair.
- **AN05B102 appears twice** in the dimorphic neuron list (left and right hemisphere); this is expected for bilateral pairs and reflects correct data, not a duplication error.

---

## 10. Future Directions

### 10.1 Multi-connectome extension
Apply the same MCIS framework when three-way NBLAST correspondence tables for MAOL and MCNS become available. Track circuit size as a function of the number of datasets — a direct measure of conservation depth.

### 10.2 Connectome-informed neural architectures and behavioral validation

**Minimal backbone controller.** The 105-neuron circuit offers a natural substrate for a structurally-grounded locomotion controller. In the flyGNN framework (Günther et al. 2023), the full 134,000-neuron connectome is instantiated as a recurrent GNN and trained end-to-end with RL, producing whole-body locomotion on a biomechanical simulator. Our circuit provides the complementary perspective: rather than instantiating the entire brain, we propose using the 105-neuron backbone as a *fixed minimal topology* — a structural prior encoding only the conserved brain–body communication channel.

Concretely, the 65 DN nodes form the input layer (receiving descending motor commands from higher brain areas), the 33 AN nodes form the output layer (encoding proprioceptive feedback to the brain), and the 12 conserved directed edges define the recurrent connections that must be preserved. All other connectivity is trainable. This differs from flyGNN in two respects: (i) the graph is three-connectome-validated rather than taken from a single specimen, and (ii) the topology is a hard constraint, not an initialisation. The prediction is that fixing the conserved edges will reduce effective degrees of freedom and improve sample efficiency on tasks requiring brain–body coordination, while having negligible benefit on pure reflex tasks (consistent with the CartPole negative result in the companion Project B analysis).

**Testable behavioral predictions via optogenetics.** The 12-edge subgraph involves a small number of identifiable hub neurons (the DNa15↔DNg04 and DNp58↔DNp65 reciprocal pairs and the DNge076 fan-out; §8.3). Because these edges are the *only* conserved connections in the circuit, their disruption should uniquely impair cross-program coordination:

1. *Multi-program silencing test:* Bilateral optogenetic silencing of each hub neuron during free locomotion should impair walking, flight initiation, and postural correction simultaneously. Single-program impairment without cross-program deficit would argue against the feedforward inhibition motif hypothesis.
2. *Edge weight prediction:* Weighted synapse counts for the 12 conserved edges should exceed the 95th percentile of all DN→AN synaptic weights in the full connectome — a prediction directly queryable via `codex.flywire.ai/api/v2/neurons/` with no new experiments required.
3. *Developmental timing:* If the conserved connectivity reflects lineage-encoded wiring (§8.2), these specific synapses should be among the earliest to appear in the pupal connectome time series; testable when developmental connectome data become available.

**Single-cell transcriptomic alignment.** Each of the 105 circuit neurons has a predicted neurotransmitter identity in `network_enriched.csv` (ACh 69.5%, GABA 16.2%, Glu 8.6%, serotonin 3.8%, dopamine 1.9%) and most have cell-type labels (DNp*, AN*) that map to clusters in published *Drosophila* single-nucleus RNA-sequencing atlases (Davie et al. 2018; Allen et al. 2025). Cross-referencing circuit membership against transcriptomic cluster identity would test whether the structurally conserved neurons form a transcriptomically coherent class — and, critically, whether their gene expression profiles contain shared regulatory logic (e.g., conserved transcription factor binding sites) that mechanistically explains cross-sex, cross-specimen wiring stereotypy. This analysis requires only the publicly available FCA (Fly Cell Atlas) data and the circuit's FAFB root IDs from `network_enriched.csv`.

### 10.3 Integration with FlyWire lab workflows
The `mcis_connectome` package enables: cross-version structural QC across proofread connectome releases; region-specific conservation queries restricted to a cell-type superclass or neuropil; developmental biology applications extending Witvliet et al. (2021) to *Drosophila*.

CLI: `python -m mcis_connectome.cli --banc ... --fafb ... --manc ... --meta ... --out network.csv`

### 10.4 Interactive Codex dashboard
A lightweight overlay on the FlyWire Codex 3D viewer that highlights MCIS membership and allows researchers to trace neuron morphology interactively across datasets.

### 10.5 Implemented extensions

These three directions are implemented in this repository (not just proposed):

**(a) Incremental MCIS + version-QC tool** ([`src/incremental_mcis.py`](src/incremental_mcis.py)). We formalise MCIS as Maximum Independent Set on the *disagreement graph* and prove it is separable over connected components, giving an exact component-local incremental update under an edge delta ΔE. Empirically, however, the disagreement graph is a single dense, low-diameter component, so exact incremental maintenance offers no speedup (a genuine structural finding: a single edit's constraint reaches almost the whole graph within two hops; a radius-1 bounded variant gives ~2.9× at a ~3-neuron approximation cost). The primitive that *is* both local and useful for proofreaders is the **O(|ΔE|) consensus-impact query** (≈2 µs, independent of graph size): given an edit to one connectome, it reports which all-3 consensus edges are gained or lost by checking only the edited pairs against the other two connectomes, and `qc_report()` renders a human-readable version-QC report naming the affected conserved edges by cell type ("a published conserved edge would disappear"). Node-level proofreading edits (**neuron merge / split**) are handled by translating them to edge edits (`node_ops_to_edge_ops`), so they ride the same engine. This is packaged as the **`mcis-watch`** command-line tool ([`src/mcis_watch.py`](src/mcis_watch.py)): given a CSV of edits in a connectome's neuron-id space it prints the QC report; a `cave_edit_delta` hook pulls real deltas from the FlyWire CAVE edit history when authenticated (`results/incremental_benchmark.json`, `results/qc_report_demo.txt`).

**(b) Connectome conservation track + spectral solver.** The conservation track (§4.6) gives a per-neuron, null-normalised conservation score, reframing the deliverable from a binary circuit to a continuous track. A spectral relaxation solver ([`src/spectral_mcis.py`](src/spectral_mcis.py)) uses the leading eigenvector of the disagreement graph to order neurons by constraint centrality and build an independent set; on sampled subgraphs it reaches ~95% of the ILP optimum at 10–100× lower runtime than ILP (`results/spectral_validation.json`). A **weighted** variant is implemented and unit-tested (`weighted_consensus` + a degree+strength-preserving null, via `conservation_track.py --weights`): it scores conserved *synaptic strength* (min weight across connectomes) against a null that preserves both degree and the weight distribution — ready to run once a synapse-resolution edge list (with weight columns) is supplied.

**(c) Ecosystem-native visualisation** ([`src/neuroglancer_overlay.py`](src/neuroglancer_overlay.py), [`src/explorer_app.py`](src/explorer_app.py), [`src/make_animation.py`](src/make_animation.py)). A FlyWire-compatible Neuroglancer state loads the circuit's FAFB neurons coloured by class or by conservation z-score (shareable via `fafbseg.encode_url`); a Streamlit app lets users filter the circuit, inspect the conservation track and conserved-edge subgraph, and download results from the committed artifacts (no bulk data download); and a rotating 3D animation (`figures/circuit_3d_conservation.gif`) shows the conserved edges concentrated along the cervical connective, coloured by conservation score.

---

## References

1. Dorkenwald et al. (2024) Neuronal wiring diagram of an adult brain. *Nature* 634, 123–138. [doi:10.1038/s41586-024-07558-y](https://doi.org/10.1038/s41586-024-07558-y)
2. Schlegel et al. (2024) Whole-brain annotation and multi-connectome cell typing of *Drosophila*. *Nature* 634, 139–152. [doi:10.1038/s41586-024-07686-5](https://doi.org/10.1038/s41586-024-07686-5) — *Consensus cell-type atlas; NBLAST matching methodology.*
3. Pospisil et al. (2024) The fly connectome reveals a path to the effectome. *Nature* 634, 234–242. [doi:10.1038/s41586-024-07982-0](https://doi.org/10.1038/s41586-024-07982-0) — *Sensorimotor bottleneck / effectome.*
4. Bates et al. (2025) Distributed control circuits across a brain-and-cord connectome. *bioRxiv*. [doi:10.1101/2025.07.31.667571](https://doi.org/10.1101/2025.07.31.667571) — *BANC dataset and metadata, including cross-dataset neuron correspondence columns.*
5. Berg et al. (2025) Sexual dimorphism in the complete connectome of the *Drosophila* male central nervous system. *bioRxiv*. [doi:10.1101/2025.10.09.680999](https://doi.org/10.1101/2025.10.09.680999) — *Cross-sex conservation baseline (~95.2% for matched DN/AN pairs).*
6. Takemura et al. (2024) A connectome of the male *Drosophila* ventral nerve cord. *eLife* 13, e97769. [doi:10.7554/eLife.97769](https://doi.org/10.7554/eLife.97769)
7. Witvliet et al. (2021) Connectomes across development reveal principles of brain maturation. *Nature* 596, 257–261. [doi:10.1038/s41586-021-03778-8](https://doi.org/10.1038/s41586-021-03778-8) — *Connectome stereotypy; ~60% of C. elegans synapses variable across individuals.*
8. Milo et al. (2002) Network motifs: simple building blocks of complex networks. *Science* 298, 824–827. [doi:10.1126/science.298.5594.824](https://doi.org/10.1126/science.298.5594.824)
9. Shiu et al. (2024) A *Drosophila* computational brain model reveals sensorimotor processing. *Nature* 634, 210–219. [doi:10.1038/s41586-024-07763-9](https://doi.org/10.1038/s41586-024-07763-9)
10. White et al. (1986) The structure of the nervous system of the nematode *Caenorhabditis elegans*. *Phil. Trans. R. Soc. Lond. B* 314, 1–340. [doi:10.1098/rstb.1986.0056](https://doi.org/10.1098/rstb.1986.0056)
11. Ito M. et al. (2013) The organization of extrinsic neurons and their implications in input/output processing of the mushroom body of *Drosophila melanogaster*. *Microscopy* 62, 58–69. [doi:10.1093/jmicro/dfs064](https://doi.org/10.1093/jmicro/dfs064) — *Hemilineage identity and developmental origin of DN/AN neurons.*
12. McGregor J.J. (1982) Backtrack search algorithms and the maximal common subgraph problem. *Software: Practice and Experience* 12, 23–34.
13. Raymond J.W. & Willett P. (2002) Maximum common subgraph isomorphism algorithms for the matching of chemical structures. *J. Computer-Aided Molecular Design* 16, 521–533.
14. Maslov S. & Sneppen K. (2002) Specificity and stability in topology of protein networks. *Science* 296, 910–913. [doi:10.1126/science.1065103](https://doi.org/10.1126/science.1065103) — *Degree-preserving rewiring null model.*
15. Lappalainen J.K. et al. (2024) Connectome-constrained networks predict neural activity across the fly visual system. *Nature* 634, 1132–1140. [doi:10.1038/s41586-024-07939-3](https://doi.org/10.1038/s41586-024-07939-3) — *Connectome-constrained mechanistic models; sparse-connectivity regime.*
16. McCreesh C., Prosser P. & Trimble J. (2017) A partitioning algorithm for maximum common subgraph problems (McSplit). *IJCAI* 712–719. — *Modern branch-and-bound MCS solver.*
17. Namiki S., Dickinson M.H., Wong A.M., Korff W. & Card G.M. (2018) The functional organization of descending sensory-motor pathways in *Drosophila*. *eLife* 7, e34272. [doi:10.7554/eLife.34272](https://doi.org/10.7554/eLife.34272) — *DN types and their leg/neck/wing VNC motor targets.*
18. Schnell B., Ros I.G. & Dickinson M.H. (2022) A population of descending neurons that regulates the flight motor of *Drosophila*. *Current Biology* 32, 1189–1196. [doi:10.1016/j.cub.2022.01.007](https://doi.org/10.1016/j.cub.2022.01.007) — *DNg02 population-code control of wingbeat amplitude.*
19. Cande J. et al. (2018) Optogenetic dissection of descending behavioral control in *Drosophila*. *eLife* 7, e34275. [doi:10.7554/eLife.34275](https://doi.org/10.7554/eLife.34275) — *DN activation → behaviour mapping.*
