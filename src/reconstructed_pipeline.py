"""
DEPRECATED — kept only as a redirect.
=====================================

This script previously contained a *second*, slightly different MCIS
implementation (it did not restrict to the giant consensus component,
removed the top-5 disagreement nodes per iteration, and had a separate
"sensorimotor-only" search that could silently overwrite network.csv).
That divergence meant the committed result and the packaged solver were not
guaranteed to agree.

The single canonical pipeline is now:

    src/run_analysis.py        (driver, multi-seed + null models, writes
                                network.csv / network_enriched.csv /
                                results/canonical_results.json)
    src/mcis_connectome/solver.py   (MCISSolver — the one algorithm)

Original method assumptions (unchanged, for reference):
  1. Use all NBLAST-matched neurons per cell type (many types have L+R pairs).
  2. Edge existence (not weight) defines circuit structure.
  3. NBLAST correspondence from BANC metadata is ground truth
     (expert annotation, Schlegel 2024 / Bates 2025 methodology).

Run instead:
    MCIS_DATA_DIR=/path/to/data python src/run_analysis.py --seeds 100
"""
import sys

if __name__ == "__main__":
    sys.stderr.write(__doc__ + "\nDelegating to src/run_analysis.py ...\n\n")
    from run_analysis import main  # noqa: E402
    main()
