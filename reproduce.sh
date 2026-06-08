#!/usr/bin/env bash
# =============================================================================
# reproduce.sh — one command to regenerate every reported number from scratch.
# =============================================================================
#
# Determinism contract
# --------------------
#   * The headline result (N, the conserved circuit, all null Z-scores) is
#     FULLY DETERMINISTIC: every stochastic step uses a fixed seed range
#     (greedy multi-start = seeds 0..99; nulls = fixed seed offsets). Re-running
#     this script on the same input files reproduces results/*.json bit-for-bit,
#     on any machine.
#   * The only non-determinism is CBC's branch-and-bound *time limit* in
#     exact_full_mis.py: the certified lower bound is stable, but the exact
#     upper bound found within the time budget can vary with CPU speed. The
#     greedy lower bound (N) does not.
#
# Inputs (point MCIS_DATA_DIR at the directory containing them):
#     banc_meta.feather   banc_626_edge_list.csv (or "... (2).csv")
#     fafb_783_edge_list.csv   manc_1.2.1_edge_list.csv
#     fafb_annotations.tsv  (enrichment only)
#
# Usage:
#     MCIS_DATA_DIR=/path/to/data ./reproduce.sh           # full pipeline
#     MCIS_DATA_DIR=/path/to/data ./reproduce.sh --quick   # skip slow certs
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

QUICK=0
[[ "${1:-}" == "--quick" ]] && QUICK=1

if [[ -z "${MCIS_DATA_DIR:-}" ]]; then
  echo "ERROR: set MCIS_DATA_DIR to the directory holding the edge lists + feather." >&2
  exit 1
fi
export MCIS_DATA_DIR
export PYTHONUNBUFFERED=1
mkdir -p results figures

step () { echo; echo "===> $*"; }

# ---- 1. Canonical result: N, the circuit, the three null models -------------
step "run_analysis.py  (N, seed distribution, correspondence + degree nulls)"
python3 src/run_analysis.py --seeds 100 --null-trials 30

# ---- 2. Derived composition / enrichment / sexual-conservation stats ---------
step "derived_stats.py  (composition, NT, descending/ascending enrichment)"
python3 src/derived_stats.py

step "confidence_tiers.py  (N vs NBLAST-confidence tier)"
python3 src/confidence_tiers.py

# ---- 3. Conservation track (now well-mixed degree null; see §4.6) -----------
step "conservation_track.py  (beyond-degree edge conservation, mixed null)"
python3 src/conservation_track.py --null-trials 50

# ---- 4. Null-model sensitivity (mixing sweep; rewire-one vs all) -------------
step "null_sensitivity.py  (swap-count convergence + rewire-one/all)"
python3 src/null_sensitivity.py --trials 15 --multipliers 0.5 1.0 3.0 10.0

# ---- 5. Heuristic validation: optimality gap under honest sampling ----------
step "exact_ilp.py  (greedy vs ILP gap; uniform sampler)"
python3 src/exact_ilp.py --sampler uniform --out ilp_validation.json
step "exact_ilp.py  (gap under HARD disagreement-ego sampling)"
python3 src/exact_ilp.py --sampler disagreement_ego --out ilp_validation_ego.json
step "exact_ilp.py  (gap under HARD degree-stratified sampling)"
python3 src/exact_ilp.py --sampler degree_stratified --out ilp_validation_stratified.json

step "worstcase_greedy.py  (synthetic systematic-underestimation demo)"
python3 src/worstcase_greedy.py

step "spectral_mcis.py  (spectral vs greedy vs ILP quality/runtime)"
python3 src/spectral_mcis.py

# ---- 6. Robustness to reconstruction error ----------------------------------
step "stringency_sweep.py  (N vs per-connectome edge perturbation)"
python3 src/stringency_sweep.py

# ---- 7. Full-graph exact MIS certificate (slow; skipped with --quick) -------
step "exact_full_mis.py  (full 987-node MIS certificate: LB <= N <= UB)"
python3 src/exact_full_mis.py --seeds 40   # clique-cover UB is instant; LB = greedy multi-start

# ---- 8. Figures + abstract ---------------------------------------------------
step "figures + abstract"
python3 src/figures_extra.py        || echo "   (figures_extra skipped)"
python3 src/regenerate_figure7.py   || echo "   (figure7 skipped)"
python3 src/make_abstract.py        || echo "   (abstract skipped)"

echo
echo "===> DONE. All numbers are in results/*.json ; figures in figures/."
echo "     See results/manifest.md for the number -> script -> output mapping."
