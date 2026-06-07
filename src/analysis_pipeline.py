"""
DEPRECATED — early cell-type-level exploration. Do NOT run.
===========================================================

This was the first attempt at the cross-connectome MCIS, using *cell-type*
annotations as correspondence keys (cell-type-level graphs across 5 datasets).
It was superseded by the individual-neuron NBLAST-correspondence approach and
is kept only for provenance. It uses a different method and historically wrote
its own `network.csv`, so it must not be executed alongside the canonical
pipeline.

Use instead:
    MCIS_DATA_DIR=/path/to/data python src/run_analysis.py --seeds 100

The original cell-type-level code is preserved in git history (see commits
prior to the integrity remediation). This shim intentionally does nothing.
"""
import sys

if __name__ == "__main__":
    sys.exit(
        "analysis_pipeline.py is deprecated and disabled to prevent it from "
        "overwriting canonical outputs.\nRun src/run_analysis.py instead."
    )
