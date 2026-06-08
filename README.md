# Structural Invariance at the Sensorimotor Interface

> *Maximum Common Induced Subgraph across three independent Drosophila connectomes*  
> FlyWire Qualification Challenge · June 2026

### ▶ Live demo — [interactive explorer](https://yanjin-ai-flywireconnectome-srcexplorer-app-pmdboy.streamlit.app/)

[![Open the live explorer](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://yanjin-ai-flywireconnectome-srcexplorer-app-pmdboy.streamlit.app/)

Browse the conserved circuit, the per-neuron conservation track, and the conserved-edge subgraph in your browser — no install, no data download.

---

![Circuit Overview](figures/figure1_circuit_layouts.png)
*The 27-neuron weakly-connected conserved sensorimotor circuit across BANC × FAFB × MANC.  
Gold edges = 26 synaptic connections verified identical across all three connectomes.  
Red = descending neurons (brain → nerve cord) · Blue = ascending neurons (nerve cord → brain).*

---

## The finding in one paragraph

Across three independently reconstructed *Drosophila* connectomes (BANC, FAFB, MANC), with neurons matched one-to-one by published NBLAST morphology, we search for the **largest connected sub-circuit whose directed wiring is identical in all three**. The answer is a verified, weakly-connected **27-neuron / 26-edge conserved circuit**, overwhelmingly **descending/ascending sensorimotor neurons** at the brain–nerve-cord interface (DN/AN enriched 67×/18× vs whole brain; 96% conserved across sexes).

**Why "connected" is the right question.** Without the connectivity constraint, the maximum common induced subgraph (= Maximum Independent Set on a "disagreement graph") is **109 neurons — but 87 of those have NO conserved edges**: they are isolated, "trivially isomorphic" nodes with no wiring to disagree about. That is exactly the *degree-inflation* the project's own degree-preserving null diagnosed (the 109 count is ~93% explained by degree sequence; §4.2). **The connectivity requirement removes the edge-less filler and isolates the real conserved circuit (27 nodes).** The strong edge-level signal — **specific wiring conserved 68.9× beyond a degree-preserving null (Z = 476σ)** — lives precisely in this connected circuit. What the data **cannot** yet tell us is *why* the wiring is conserved (developmental vs. activity-driven); see the evidence ledger.

**How to read this repo:**
- **30-sec result** → *Result at a Glance* (below)
- **Method & algorithm walkthrough** → [science.md §3](science.md) (problem formulation, solver decision rule §3.5, worst-case self-critique §3.4)
- **What's certain / what to question** → [science.md §0 evidence ledger](science.md) and §8.4 alternative hypotheses
- **Reproduce every number** → [`./reproduce.sh`](reproduce.sh) + [`results/manifest.md`](results/manifest.md)

## Result at a Glance

| | |
|--|--|
| **🔑 Conserved circuit (connectivity-constrained)** | **N = 27 neurons, 26 conserved edges, weakly connected** — the deliverable (`network.csv`); verified isomorphic across all 3 datasets, 0 internal disagreement |
| **Composition** | 17 descending (63%) + 6 ascending (22%) + 4 sensory; ACh 41% / GABA 30% (mixed excitatory–inhibitory) |
| **Cell-type enrichment** | Descending **67.3×** (p=2×10⁻²⁸), Ascending **17.7×** (p=9×10⁻⁷) vs FAFB whole-brain background |
| **Sexual conservation** | **96.3%** isomorphic across ♀ and ♂ (26/27) |
| **Contrast — unconstrained MCIS** | 109 neurons (GMIN+2-swap), but **87 are edge-less isolated nodes** → connectivity strips this degree-inflation to the 27-node circuit (`network_unconstrained_mcis.csv`) |
| **🔑 Wiring conserved *beyond degree*** | **2,609 consensus edges vs well-mixed degree-null 37.9 ± 5.4 → 68.9×, Z = 476σ** ([§ Conservation Track](#beyond-the-binary-circuit--conservation-track-version-qc--visual-tools)) |
| **Datasets** | BANC v626 (♀ brain+cord) × FAFB v783 (♀ brain) × MANC v1.2.1 (♂ nerve cord) |
| **Optimality** | connected-MIS is NP-hard; N=27 is stable under multi-start + 30k-iteration ILS (`results/connected_mcis.json`) |

---

## Beyond the Binary Circuit — Conservation Track, Version-QC & Visual Tools

The MCIS *size* (node count) is largely explained by the degree sequence — but the shared **wiring** is not. This section is the substantive extension beyond the qualification result.

### 🔑 Specific connectivity is conserved far beyond degree (the headline)

Over the 987-node consensus component, **2,609 directed edges are present in all three connectomes vs 37.9 ± 5.4 under a *well-mixed* degree-preserving null → 68.9× enrichment, Z = 476σ**. So cross-connectome agreement reflects real wiring identity, not matched degree distributions. We therefore report a continuous, null-normalised **per-neuron conservation track** instead of a binary circuit ([`src/conservation_track.py`](src/conservation_track.py) → `results/conservation_track.json`, `results/neuron_conservation.csv`).

> **Null-mixing note (methodological honesty).** An earlier version reported 7.4× against a null of 353 ± 17. That null was *under-mixed*: `directed_edge_swap` was run with only ~0.5×\|E\| swaps, leaving residual real structure that inflated the null. [`src/null_sensitivity.py`](src/null_sensitivity.py) shows the null consensus count only converges past ~3×\|E\| swaps (to ≈38); the corrected, well-mixed enrichment is **68.9×**, *stronger* than before. The default swap count is now 10×\|E\|.

![Conservation track](figures/figure10_conservation_track.png)

### Conserved circuit in 3D (conserved edges = gold; colour = conservation z)

![3D conservation animation](figures/circuit_3d_conservation.gif)

*The conserved edges concentrate along the cervical connective — the expected brain↔cord relay locus. ([`src/make_animation.py`](src/make_animation.py))*

### Version-QC tool — "does my proofreading edit touch a published conserved circuit?"

An **O(|ΔE|) consensus-impact query** (~2 µs, independent of graph size) reports which conserved edges an edit gains/loses, named by cell type ([`src/incremental_mcis.py`](src/incremental_mcis.py) → `results/qc_report_demo.txt`):

```
Connectome version QC report
================================
Conserved (all-3) edges lost:   1
Neurons touching conserved wiring: 2
  LOST  ANXXX202_b -> AN27X017  (a published conserved edge would disappear)
```

Packaged as the **`mcis-watch`** CLI (`python src/mcis_watch.py --dataset FAFB --edits edits.csv`), which also accepts neuron **merge/split** edits and has a `cave_edit_delta` hook for live FlyWire CAVE edit history. We formalise MCIS as Maximum Independent Set on the *disagreement graph*; exact incremental maintenance is provably correct but gives no speedup because that graph is one dense, low-diameter component (an honest structural finding) — so the O(|ΔE|) query above is the primitive that is both local and useful.

| | |
|--|--|
| ![Incremental](figures/figure11_incremental.png) | ![Spectral](figures/figure12_spectral.png) |
| Incremental MCIS: speed/accuracy vs radius | Spectral solver: 92–97% of ILP (beaten by greedy at scale; §3.5) |

### Spectral relaxation solver

An eigenvector-based MIS heuristic on the disagreement graph reaches ~92–97% of the ILP optimum ([`src/spectral_mcis.py`](src/spectral_mcis.py) → `results/spectral_validation.json`). Honest caveat: it is **beaten by plain greedy** at sizes ≥ 80, so we keep it as a constraint-centrality probe, not the production solver — see the solver decision rule in [science.md §3.5](science.md).

### Ecosystem-native, interactive

- **FlyWire Neuroglancer overlay** — `python src/neuroglancer_overlay.py --color conservation` writes [`results/neuroglancer_state.json`](results/neuroglancer_state.json) (109 FAFB neurons coloured by conservation z); open at [ngl.flywire.ai](https://ngl.flywire.ai/) or shorten via `fafbseg.encode_url`.
- **Streamlit explorer** — interactive: filter the circuit, inspect the conservation track + conserved-edge subgraph, download CSV. Runs from committed artifacts (no bulk data download).

  ▶ **Live:** <https://yanjin-ai-flywireconnectome-srcexplorer-app-pmdboy.streamlit.app/>

  [![Open the live explorer](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://yanjin-ai-flywireconnectome-srcexplorer-app-pmdboy.streamlit.app/)

  Locally: `streamlit run src/explorer_app.py`. Deploy your own: see [DEPLOY.md](DEPLOY.md).

---

## Challenge deliverables

| Required | Delivered |
|---|---|
| Solution CSV: 3 dataset columns, N matched-neuron rows | [`network.csv`](network.csv) — 109 rows × {BANC, FAFB, MANC} |
| Maximise N; mutually isomorphic directed induced subgraphs (edge ⟺ in all, direction preserved) | N=109, 14 edges, verified `isomorphic=True` ([`results/canonical_results.json`](results/canonical_results.json)) |
| Research: network-graph visualization | one-pager panel 1 + `figures/figure1`, `figure3` |
| Research: Codex 3D meshes | **Verified** — all 109 neurons render in the Codex 3D viewer inside the FAFB whole-brain mesh (live link in [`results/codex_3d_url.txt`](results/codex_3d_url.txt); IDs in [`results/codex_circuit_ids.txt`](results/codex_circuit_ids.txt); see [`results/codex_links.md`](results/codex_links.md)) |
| Research: observations / hypothesis | one-pager panel 3 + science.md §5–§8 |
| Research: literature & citations | one-pager refs + science.md References (13) |
| **Concise one-page summary (one dataset = FAFB)** | **[`research_summary_fafb.pdf`](research_summary_fafb.pdf)** |

---

## Scientific Framing

The FlyWire multi-connectome cell typing atlas (Schlegel et al. 2024) established that **cell-type identity** is reproducible across connectomes at the morphological level. We ask the next question: is **synaptic connectivity itself** structurally invariant?

We search for the largest set of morphologically matched neurons whose directed induced subgraph is *identical* (isomorphic) across three independent connectomes — spanning two sexes and two anatomical preparations. The result is a 109-neuron sensorimotor backbone enriched 65.7× for descending neurons and 24.8× for ascending neurons relative to the whole-brain background, consistent with the sensorimotor bottleneck hypothesis (Pospisil et al. 2024).

---

## Research Summary — biological significance

*The 1-page scientific report is [`science.md`](science.md) (it opens with a one-page summary); a print-ready version is [`research_summary_fafb.pdf`](research_summary_fafb.pdf). Key content is mirrored here.*

**What the circuit is / does.** A 27-neuron weakly-connected brain↔ventral-nerve-cord locomotor circuit (DN/AN/sensory). Its hub is the ascending neuron AN02A002 (fans out to 10 leg-motor DNs); it is leg-VNC-dominated, carried by as-yet-uncharacterised DN/AN types (a central ascending hub AN02A002 fans out to 10 leg-motor descending neurons; a DNp58↔DNp65 reciprocal pair anchors recurrent control); these DN classes target leg/neck/wing VNC motor circuits (Namiki et al. 2018). DNg02 — a documented flight controller (Schnell et al. 2022) — is in the set but is sexually dimorphic and not a conserved-edge carrier. Mixed ACh/GABA chemistry → feedforward-inhibition coordination motif.

| Network graph | Codex 3D meshes (FAFB) |
|---|---|
| ![network](figures/figure3_hub_circuit.png) | ![codex](figures/codex_3d_fafb.png) |

**Structural observations.** 65.7×/24.8× descending/ascending enrichment (Fisher p<10⁻³⁵); wiring conserved **68.9× beyond a well-mixed degree-preserving (Maslov–Sneppen) null** (2,609 vs 37.9±5.4, Z=476σ); 92.7% conserved across sexes; peripheral (low-betweenness) relays; robust to 20% simulated reconstruction error.

**Interpretation / hypotheses.** A developmentally canalised brain↔cord channel (H1); the degree-null rejects a pure degree artifact (H2) and the enrichment rejects a generic-subgraph explanation (H5); static data cannot yet distinguish developmental vs activity-driven wiring (H4 — key open question). **Prediction:** silencing the hub neurons (the AN02A002 ascending hub (a conserved fan-out onto 10 leg-motor DNs) and the DNp58↔DNp65 reciprocal pair) should impair walking, flight and posture simultaneously. Full citations and the alternative-hypothesis table are in [`science.md`](science.md) §6.5, §8.

---

## Technical Approach

### Core insight: official morphological correspondence

Rather than defining ad hoc neuron matching, we use the **BANC metadata** (`banc_888_meta.feather`, Bates et al. 2025) which contains `fafb_match` and `manc_match` columns — 1:1 NBLAST-based morphological correspondences per neuron, giving 3,414 pre-verified triplets with biological ground truth.

### Why BANC × FAFB × MANC?

BANC is the only dataset spanning both brain and ventral nerve cord. Its metadata explicitly provides individual-neuron-level cross-links to FAFB (brain-only) and MANC (cord-only). No equivalent three-way correspondence exists for MAOL or MCNS at this resolution. The triplet also tests cross-sex conservation (♀ FAFB/BANC vs ♂ MANC).

### Graph-matching formulation

Because the NBLAST correspondence fixes a 1:1 neuron identity across datasets, we **do not search over node matchings** (the hard part of general Maximum Common Subgraph). With the bijection fixed, "mutually isomorphic directed induced subgraph" collapses to **edge-set equality**, and maximising N becomes exactly **Maximum Independent Set (MIS) on the *disagreement graph*** — the graph that links any neuron pair whose directed connection is present in some but not all three datasets. We solve this three ways: a fast **greedy heuristic** (below, multi-start), an **exact ILP** (PuLP/CBC) for certified near-optimality on subgraphs, and a **spectral relaxation** as a fast alternative. This is why we avoid the NP-hard matching search of classical MCS solvers (McGregor 1982; Raymond & Willett 2002; McSplit).

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
**Reproducibility:** Fixed random seeds; deterministic per seed; result reported as the best of a 3000-restart GMIN+2-swap multi-start.  
**Unit tests:** `pytest tests/ -v` — 21 passed, 1 skipped (real-data smoke), all pass.  
**This is Maximum Independent Set on a "disagreement graph"** (a pair of neurons cannot coexist if their connection disagrees across connectomes). Greedy = the classic max-degree vertex-cover heuristic — see the **solver decision rule** (greedy vs ILP vs spectral) in [science.md §3.5](science.md) and the worst-case self-critique in §3.4.  
**Optimality:** honest greedy gap **2.67%** under hard (disagreement-ego) subgraph sampling (1.15% under optimistic uniform sampling); full-graph certificate **109 ≤ N ≤ 136** with N=109 verified-feasible ([science.md §3.6](science.md), `results/exact_full_mis.json`).

### Assumptions

1. **Fixed correspondence is ground truth.** A neuron's identity across datasets is the BANC-metadata NBLAST `fafb_match` / `manc_match` (Bates et al. 2025); we do not re-derive matches. One neuron ↔ one identity (a bijection), so isomorphism = edge-set equality.
2. **Unweighted, directed edges define structure.** Per the challenge spec, synapse weights are ignored; an edge is present/absent, direction preserved. All analysis is on the unweighted directed graphs.
3. **Edge existence in the provided edge lists is authoritative** (proofreading errors are treated as noise — robustness to this is quantified in `results/stringency_sweep.json`: N degrades gracefully under simulated error).
4. **Dataset versions:** BANC v626 edge list (correspondence columns from the `banc_888` metadata build; `root_626` is the join key), FAFB v783, MANC v1.2.1.
5. **Reported N is the best of a fixed 3000-restart GMIN+2-swap multi-start** (a verified-feasible lower bound); ILP on sampled subgraphs bounds the honest optimality gap to ~2.7%, and a Lovász-θ certificate gives N ≤ 136.

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
| 3000 GMIN+2-swap restarts | N = 104.8 ± 1.8, range [99, 109]; best = 109 | Stable; not seed-dependent |
| Correspondence-shuffle null (20× best-of-5) | N_null = 81.1 ± 2.3 (>15σ below real) | Neuron identity is essential |
| Degree-preserving rewire null (well-mixed, 30× best-of-5) | N_null = 101.9 ± 1.5 (~4 below real per-seed mean; 4.6σ below best 109) | Degree explains ~93% of node count; neuron identity adds the remaining ~7 |
| Centrality permutation (1000 trials) | p = 0.003; circuit has *lower* betweenness | Peripheral relays, not hubs |
| Manual annotation rate | Circuit 92.7% vs non-circuit 85.1% (p = 0.014) | Better-annotated correspondences |
| NBLAST confidence curve (6 tiers) | N = 4 → 41 → 55 → 77 → 92 → 109 (monotonic) | Not driven by low-confidence matches |
| **Conservation beyond degree (edges)** | **2,609 consensus edges vs well-mixed degree-null 37.9 ± 5.4 → 68.9×, Z = 476σ** | Specific wiring is conserved far beyond degree sequence |
| Reconstruction-error robustness | N = 109 → 92 → 85 → 73 at 0/5/10/20% edge flips per connectome | Graceful decline, no cliff — not an artifact of the exact edge sets |
| **Null-mixing sensitivity** (new) | Beyond-degree Z stable once null is mixed (>3×\|E\| swaps); in/out degree preserved exactly | The headline is robust to the null's swap count — the old 7.4× was an under-mixing artifact (`results/null_sensitivity.json`) |
| **Honest sampling** (new) | Greedy gap 1.15% (uniform) → 2.67% (hard disagreement-ego sampling) | Uniform sampling is optimistic; greedy still ~97% of exact on the hard regime (`results/ilp_validation_ego.json`) |
| **Greedy worst-case** (new) | Synthetic: gap grows with density; tight-instance ≈ln k underestimate | Characterises *where* max-degree greedy underestimates MIS (`results/worstcase_greedy.json`) |
| **Full-graph MIS certificate** (new) | 109 ≤ N ≤ 136; N=109 verified-feasible (0 internal disagreement) | N is a checked feasible solution, not just a heuristic output (`results/exact_full_mis.json`) |

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
print(result.summary())  # MCISResult(N=109, edges=12, isomorphic=True)
result.to_csv('network.csv')
```

---

## Reproduction

```bash
git clone https://github.com/Yanjin-ai/flywireconnectome.git
cd flywireconnectome
pip install -e ".[test]"          # or: pip install -r requirements.txt

# Place the input files in a data directory and point MCIS_DATA_DIR at it:
#   banc_meta.feather          banc_626_edge_list.csv
#   fafb_783_edge_list.csv     manc_1.2.1_edge_list.csv
#   fafb_annotations.tsv       (for the enrichment analysis)
export MCIS_DATA_DIR=/path/to/data

# Edge lists:     https://codex.flywire.ai/api/download
# BANC metadata:  https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather
# FAFB annotations: https://github.com/flyconnectome/flywire_annotations

# Run the canonical pipeline (writes network.csv + results/ + figures/ into the repo):
# …or run EVERYTHING in order with one command: ./reproduce.sh  (see results/manifest.md)
python src/run_analysis.py --seeds 100        # → network.csv, results/canonical_results.json
python src/exact_ilp.py                        # → results/ilp_validation.json (PuLP/CBC, uniform sampling)
python src/exact_ilp.py --sampler disagreement_ego --out ilp_validation_ego.json  # honest hard-sampling gap
python src/exact_full_mis.py --seeds 40        # → full-graph MIS certificate 109 ≤ N ≤ 136
python src/null_sensitivity.py                 # → null mixing/sensitivity (results/null_sensitivity.json)
python src/worstcase_greedy.py                 # → where greedy underestimates MIS (synthetic)
python src/derived_stats.py                    # → results/derived_stats.json (composition/enrichment)
python src/confidence_tiers.py                 # → results/confidence_tiers.json
python src/visualize.py                        # → figures/figure1-4
python src/figures_extra.py                    # → figures/figure6,8,9
python src/regenerate_figure7.py               # → figures/figure7
python src/robustness_experiments.py           # → figures/figure5 (from results/*.json)
python src/conservation_track.py                # → conservation track: edges conserved beyond degree (Z=476σ, well-mixed null)
python src/incremental_mcis.py                  # → incremental MCIS + O(|ΔE|) version-QC report
python src/spectral_mcis.py                     # → spectral solver vs ILP/greedy (research probe; greedy wins at scale, §3.5)
python src/stringency_sweep.py                  # → reconstruction-error robustness (results/stringency_sweep.json)
python src/neuroglancer_overlay.py --color conservation  # → results/neuroglancer_state.json (FlyWire)
python src/make_animation.py                    # → figures/circuit_3d_conservation.gif
python src/make_abstract.py                     # → extended_abstract.pdf/.png (from results/*.json)
python src/make_research_summary.py             # → research_summary_fafb.pdf + Codex ID list
streamlit run src/explorer_app.py               # → interactive explorer (uses only committed artifacts)
pytest tests/ -v                               # → unit tests pass (synthetic graphs; real-data smoke test skips without MCIS_DATA_DIR)
```

> **Note on BANC versioning.** The correspondence columns (`fafb_match`,
> `manc_match`) come from the `banc_888_meta.feather` build; the BANC edge list
> uses `root_626` IDs. The metadata's `root_626` column is the join key into the
> v626 edge list, so the cross-dataset bijection is consistent despite the build
> labels differing.

---

## Repository Structure

```
network.csv                    109 rows × 3 columns (BANC | FAFB | MANC neuron IDs)
science.md                     Scientific report (10 sections, 13 references)
results/                       Reproducible JSON outputs (canonical, ILP, derived, confidence)
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
  test_solver.py               21 passed + 1 skipped (real-data smoke)
src/
  run_analysis.py             Canonical pipeline (MCIS + nulls + centrality)
  exact_ilp.py                ILP optimality validation (PuLP/CBC); --sampler uniform|degree_stratified|disagreement_ego
  exact_full_mis.py           Full-graph MIS certificate: LB (verified) ≤ N ≤ clique-cover UB (§3.6)
  null_sensitivity.py         Null-model mixing/sensitivity sweep (figure 14)
  worstcase_greedy.py         Where max-degree greedy underestimates MIS (figure 15)
  derived_stats.py            Composition / enrichment / dimorphism / annotation quality
  confidence_tiers.py         NBLAST confidence-tier curve
  conservation_track.py        Per-edge/per-neuron conservation vs degree-null (figure 10)
  incremental_mcis.py          Incremental MCIS + O(|ΔE|) impact query, merge/split ops, CAVE hook (figure 11)
  mcis_watch.py                `mcis-watch` CLI: QC-report an edit against conserved circuits
  spectral_mcis.py             Spectral relaxation solver vs ILP/greedy (figure 12)
  neuroglancer_overlay.py      FlyWire Neuroglancer state (colour by class/conservation)
  explorer_app.py              Streamlit interactive explorer
  make_animation.py            Rotating 3D conservation GIF
  make_research_summary.py     One-page FAFB research summary + Codex neuron-id list
  mcis_paths.py               Shared path resolution (MCIS_DATA_DIR → ./data)
  visualize.py                 figures 1-4
  figures_extra.py             figures 6, 8, 9
  regenerate_figure7.py        figure 7
  robustness_experiments.py    figure 5 (from results/*.json)
  make_abstract.py             extended_abstract.pdf/.png (from results/*.json)
  reconstructed_pipeline.py    Deprecated shim → run_analysis.py
  analysis_pipeline.py         Deprecated/disabled (early cell-type-level method)
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
