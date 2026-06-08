"""
Per-neuron NBLAST match confidence for the 105-neuron circuit (closes the
"correspondence is a hidden single point of failure" gap).
=========================================================================

Every result rests on the BANC-metadata 1:1 neuron correspondence (NBLAST
top-1). A wrong match injects a spurious disagreement edge and can move both N
and circuit membership. A reviewer therefore asks: *how certain are the 105
matches themselves?* This script answers it with the data we actually have.

Confidence signals available in banc_meta.feather (top-1 match IDs only — no
continuous NBLAST scores or top-2, which would need the bancr R package + ~10 GB
skeletons; stated as a limitation):

  conf2  = (fafb_match == fafb_nblast_match) + (manc_match == manc_nblast_match)
           in {0,1,2}: agreement between the curated match and the automated
           NBLAST top-1, for the two datasets that define the triplet (the §4.5
           proxy, here computed PER circuit neuron, not just in aggregate).

  multi  = additionally counts agreement against every OTHER connectome's NBLAST
           top-1 that has a curated match (hemibrain / fanc / malecns): an
           independent cross-dataset corroboration of the same neuron identity.

We report the distribution over the 105 circuit neurons, compare it to the
987-node consensus background (is the circuit built on BETTER-or-WORSE matches?),
and write a per-neuron confidence column so the 105 can be audited individually.

Run:  MCIS_DATA_DIR=/path/to/data python src/match_confidence.py
Outputs: results/match_confidence.json, results/circuit_match_confidence.csv,
         figures/figure16_match_confidence.png
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(__file__))
from mcis_paths import data_dir, results_dir, figures_dir  # noqa
from run_analysis import build_solver, seed_distribution  # noqa

NBLAST_PAIRS = [
    ("fafb_match", "fafb_nblast_match"),
    ("manc_match", "manc_nblast_match"),
    ("hemibrain_match", "hemibrain_nblast_match"),
    ("fanc_match", "fanc_nblast_match"),
    ("malecns_match", "malecns_nblast_match"),
]


def agreement(row, col, nbl):
    """1 if a curated match AND an NBLAST top-1 both exist and agree; else 0.
    Returns (counts_as_available, agrees)."""
    a, b = row.get(col), row.get(nbl)
    if pd.isna(a) or pd.isna(b) or str(a) in ("", "nan") or str(b) in ("", "nan"):
        return 0, 0
    return 1, int(str(a) == str(b))


def score_rows(df):
    """conf2 (fafb+manc agreement, 0..2) and a multi-dataset corroboration
    (agreements / datasets-with-a-curated-match) for each row."""
    conf2, multi_frac, multi_n = [], [], []
    for _, r in df.iterrows():
        c2 = 0
        for col, nbl in NBLAST_PAIRS[:2]:
            _, ag = agreement(r, col, nbl)
            c2 += ag
        avail = agree = 0
        for col, nbl in NBLAST_PAIRS:
            a, ag = agreement(r, col, nbl)
            avail += a
            agree += ag
        conf2.append(c2)
        multi_n.append(avail)
        multi_frac.append(agree / avail if avail else np.nan)
    return np.array(conf2), np.array(multi_frac), np.array(multi_n)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--seeds", type=int, default=40)
    args = ap.parse_args()

    dd = Path(args.data_dir) if args.data_dir else Path(data_dir())
    solver = build_solver(dd, n_seeds=1)
    solver._load()
    trip = solver._triples

    # 105-neuron circuit (best of multi-start) as giant-local indices -> rows
    _, best = seed_distribution(solver, args.seeds)
    circ_rows = trip.iloc[[solver._gl[i] for i in sorted(best)]].reset_index(drop=True)
    # 987-node consensus background
    bg_rows = trip.iloc[[solver._gl[i] for i in range(solver._ng)]].reset_index(drop=True)

    c2_circ, mf_circ, mn_circ = score_rows(circ_rows)
    c2_bg, mf_bg, mn_bg = score_rows(bg_rows)

    def dist(c2):
        return {str(k): int((c2 == k).sum()) for k in (0, 1, 2)}

    circ_dist, bg_dist = dist(c2_circ), dist(c2_bg)
    n_circ = len(c2_circ)
    # high-confidence = both NBLAST top-1 agree (conf2 == 2)
    hi_circ = int((c2_circ == 2).sum())
    hi_bg = int((c2_bg == 2).sum())
    # enrichment of high-confidence in circuit vs background (Fisher)
    table = [[hi_circ, n_circ - hi_circ],
             [hi_bg - hi_circ, (len(c2_bg) - hi_bg) - (n_circ - hi_circ)]]
    _, p = fisher_exact(table)

    print(f"  Circuit (N={n_circ}) conf2 distribution (both/one/neither NBLAST "
          f"top-1 agree): {circ_dist}")
    print(f"  Background (N={len(c2_bg)}) conf2 distribution: {bg_dist}")
    print(f"  High-confidence (both agree): circuit {100*hi_circ/n_circ:.0f}% vs "
          f"background {100*hi_bg/len(c2_bg):.0f}%  (Fisher p = {p:.3g})")
    print(f"  Multi-dataset corroboration (mean agreement fraction over all "
          f"connectomes with a curated match): circuit "
          f"{np.nanmean(mf_circ):.2f} vs background {np.nanmean(mf_bg):.2f}")

    # per-neuron audit table
    audit = pd.DataFrame({
        "BANC": circ_rows["root_626"].astype(str),
        "cell_type": circ_rows.get("cell_type"),
        "conf2_fafb_manc_agree": c2_circ,
        "multi_dataset_agree_frac": np.round(mf_circ, 3),
        "n_datasets_corroborating": mn_circ,
    }).sort_values("conf2_fafb_manc_agree", ascending=False)
    os.makedirs(results_dir(), exist_ok=True)
    audit.to_csv(results_dir() + "circuit_match_confidence.csv", index=False)

    out = {
        "circuit_n": n_circ,
        "conf2_definition": "(fafb_match==fafb_nblast_top1)+(manc_match==manc_nblast_top1), 0..2",
        "circuit_conf2_distribution": circ_dist,
        "background_conf2_distribution": bg_dist,
        "high_conf_pct_circuit": round(100 * hi_circ / n_circ, 1),
        "high_conf_pct_background": round(100 * hi_bg / len(c2_bg), 1),
        "fisher_p_high_conf_enrichment": float(p),
        "multi_dataset_mean_agree_circuit": float(np.nanmean(mf_circ)),
        "multi_dataset_mean_agree_background": float(np.nanmean(mf_bg)),
        "caveat": ("banc_meta.feather stores top-1 match IDs only; continuous "
                   "NBLAST scores / top-1-vs-top-2 gaps need the bancr R package "
                   "+ ~10 GB skeletons. Agreement against independent connectomes' "
                   "NBLAST top-1 is the strongest available per-neuron proxy."),
    }
    with open(results_dir() + "match_confidence.json", "w") as f:
        json.dump(out, f, indent=2)

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor("white")
    labels = ["neither", "one", "both"]
    cw = [circ_dist[k] / n_circ * 100 for k in ("0", "1", "2")]
    bw = [bg_dist[k] / len(c2_bg) * 100 for k in ("0", "1", "2")]
    x = np.arange(3)
    axes[0].bar(x - 0.2, cw, 0.4, label=f"circuit (N={n_circ})", color="#d4a017")
    axes[0].bar(x + 0.2, bw, 0.4, label="987 background", color="#999")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"{l}\nNBLAST top-1 agree" for l in labels])
    axes[0].set_ylabel("% of neurons")
    axes[0].set_title("A  Per-neuron match confidence (FAFB+MANC)",
                      fontweight="bold")
    axes[0].legend()
    axes[1].hist(mf_circ[~np.isnan(mf_circ)], bins=np.linspace(0, 1, 11),
                 color="#d4a017", alpha=0.85)
    axes[1].set_xlabel("fraction of connectomes whose NBLAST top-1 agrees")
    axes[1].set_ylabel("circuit neurons")
    axes[1].set_title("B  Multi-dataset corroboration of the 105 matches",
                      fontweight="bold")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("NBLAST match confidence for the conserved circuit",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    os.makedirs(figures_dir(), exist_ok=True)
    fig.savefig(figures_dir() + "figure16_match_confidence.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print("  Wrote results/match_confidence.json, circuit_match_confidence.csv, "
          "figures/figure16_match_confidence.png")


if __name__ == "__main__":
    main()
