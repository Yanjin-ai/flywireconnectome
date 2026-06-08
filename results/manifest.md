# Results manifest — every reported number → producing script → output file

This is the provenance table for the project. **No number appears in `README.md`
or `science.md` that is not produced by one of these committed scripts writing to
`results/`.** Re-run everything with `./reproduce.sh` (see the determinism
contract at the top of that script).

| Reported quantity | Producing script | Output file (key) |
|---|---|---|
| N (best of 100-seed multi-start), conserved edges, isomorphism | `src/run_analysis.py` | `canonical_results.json` → `N_reported`, `n_conserved_edges` |
| Seed distribution (mean ± sd, range) | `src/run_analysis.py` | `canonical_results.json` → `seed_distribution` |
| Correspondence-shuffle null | `src/run_analysis.py` | `canonical_results.json` → `correspondence_shuffle_null` |
| Degree-preserving null (node-count) | `src/run_analysis.py` | `canonical_results.json` → `degree_preserving_null` |
| Centrality (betweenness, permutation p) | `src/run_analysis.py` | `canonical_results.json` → `centrality` |
| N bounds waterfall (triplets → giant → N) | `src/run_analysis.py` | `canonical_results.json` → `n_bounds` |
| Composition, neurotransmitters, sexual conservation | `src/derived_stats.py` | `derived_stats.json` |
| Descending/ascending fold-enrichment + Fisher p | `src/derived_stats.py` | `derived_stats.json` → `enrichment` |
| N vs NBLAST-confidence tier | `src/confidence_tiers.py` | `confidence_tiers.json` |
| **Beyond-degree edge conservation (mixed null)** | `src/conservation_track.py` | `conservation_track.json` → `degree_null` (enrichment, z) |
| Per-neuron conservation z-track | `src/conservation_track.py` | `neuron_conservation.csv` |
| **Null-model sensitivity** (mixing sweep, rewire-one vs all, reciprocity) | `src/null_sensitivity.py` | `null_sensitivity.json` |
| Heuristic optimality gap — uniform sampling | `src/exact_ilp.py --sampler uniform` | `ilp_validation.json` |
| Heuristic optimality gap — hard (ego) sampling | `src/exact_ilp.py --sampler disagreement_ego` | `ilp_validation_ego.json` |
| Heuristic optimality gap — hard (degree-stratified) sampling | `src/exact_ilp.py --sampler degree_stratified` | `ilp_validation_stratified.json` |
| Greedy systematic-underestimation demo | `src/worstcase_greedy.py` | `worstcase_greedy.json` |
| **Full-graph exact MIS certificate (LB ≤ N ≤ UB)** | `src/exact_full_mis.py --theta` | `exact_full_mis.json` → `certificate` |
| **Production solver beats baseline greedy (105→109)** | `src/improve_mis.py` | `improve_mis.json` |
| **Per-neuron NBLAST match confidence of the circuit** | `src/match_confidence.py` | `match_confidence.json`, `circuit_match_confidence.csv` |
| Spectral vs greedy vs ILP quality/runtime | `src/spectral_mcis.py` | `spectral_validation.json` |
| Reconstruction-error robustness (N vs edge flip p) | `src/stringency_sweep.py` | `stringency_sweep.json` |
| Incremental MCIS / consensus-impact QC benchmark | `src/incremental_mcis.py` | `incremental_benchmark.json` |
| Neuroglancer overlay state | `src/neuroglancer_overlay.py` | `neuroglancer_state.json` |

## Determinism

- **Deterministic, machine-independent:** N, the circuit (`network.csv`), all
  null means/SDs/Z-scores, every figure. Fixed seeds throughout.
- **Mildly machine-dependent:** the *upper* bound in `exact_full_mis.json` (CBC
  branch-and-bound within a wall-clock limit). The certified *lower* bound (= N)
  and the greedy result do not depend on CPU speed.
