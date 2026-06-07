"""
Robustness figure (figure5_robustness.png).
============================================

Plots the robustness panel DIRECTLY from the reproducible result files written
by the canonical pipeline, so the figure can never drift from the numbers in
science.md:

    results/canonical_results.json   (run_analysis.py)
    results/confidence_tiers.json    (confidence_tiers.py)
    results/ilp_validation.json      (exact_ilp.py)

Annotation quality (manual_cluster) is recomputed from the BANC metadata
(metadata-only, no edge-list loading required).

Run:  MCIS_DATA_DIR=/path/to/data python src/robustness_experiments.py
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mcis_paths import data_dir, repo_root
REPO = repo_root()
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(REPO, "figures", "figure5_robustness.png")
DATA_DIR = data_dir()


def load(name):
    with open(os.path.join(RESULTS, name)) as f:
        return json.load(f)


def annotation_quality():
    m = pd.read_feather(DATA_DIR + "banc_meta.feather")
    m["root_626"] = m["root_626"].astype(str)
    m = m.drop_duplicates("root_626")
    pool = m[m["fafb_match"].notna() & m["manc_match"].notna()]
    circ = set(pd.read_csv(os.path.join(REPO, "network.csv"), dtype=str)["BANC"])
    ic = pool["root_626"].isin(circ)
    nn = pool["manual_cluster"].notna()
    return 100 * nn[ic].mean(), 100 * nn[~ic].mean()


def main():
    can = load("canonical_results.json")
    tiers = load("confidence_tiers.json")["tiers"]
    ilp = load("ilp_validation.json")

    sizes = np.array(can["seed_distribution"]["sizes"])
    cs = can["correspondence_shuffle_null"]
    dp = can["degree_preserving_null"]
    cent = can["centrality"]
    best_n = can["N_reported"]

    fig = plt.figure(figsize=(20, 12))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)
    BLUE, GOLD, RED, GREEN, PURP = "#1f77b4", "#d4a017", "#d62728", "#2ca02c", "#9467bd"

    def style(ax, title):
        ax.set_title(title, fontsize=12, pad=8, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

    # A: seed distribution
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(sizes, bins=range(sizes.min(), sizes.max() + 2), color=BLUE, alpha=0.85,
            align="left")
    ax.axvline(sizes.mean(), color=GOLD, ls="--", lw=2, label=f"mean={sizes.mean():.1f}")
    ax.axvline(best_n, color=RED, lw=2, label=f"best={best_n}")
    ax.set_xlabel("MCIS size N"); ax.set_ylabel("Count"); ax.legend()
    style(ax, f"A  100-seed distribution ({sizes.mean():.1f} ± {sizes.std():.1f})")

    # B: null comparison
    ax = fig.add_subplot(gs[0, 1])
    data = [sizes, np.array(dp["sizes"]), np.array(cs["sizes"])]
    parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True)
    for pc, c in zip(parts["bodies"], [BLUE, GREEN, RED]):
        pc.set_facecolor(c); pc.set_alpha(0.6)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["Real", "Degree-\npreserving", "Correspond.\nshuffle"])
    ax.set_ylabel("MCIS size N")
    style(ax, "B  Null models (best-of-5 per trial)")
    ax.text(0.5, 0.05, f"shuffle: {cs['mean']:.1f}±{cs['std']:.1f}   "
            f"degree: {dp['mean']:.1f}±{dp['std']:.1f}",
            transform=ax.transAxes, ha="center", fontsize=9)

    # C: bound waterfall
    ax = fig.add_subplot(gs[0, 2])
    nb = can["n_bounds"]
    labels = ["matched\ntriplets", "in all 3\nedge lists", "giant\ncomponent", "MCIS N"]
    vals = [nb["matched_triplets"], nb["in_all3_edge_lists"], nb["giant_component"], best_n]
    ax.bar(labels, vals, color=[BLUE, BLUE, BLUE, GOLD], alpha=0.85)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("count"); ax.set_yscale("log")
    style(ax, "C  N bound waterfall")

    # D: ILP optimality gap
    ax = fig.add_subplot(gs[1, 0])
    szs = [t["size"] for t in ilp["tiers"]]
    gaps = [t["gap_pct"] for t in ilp["tiers"]]
    ax.plot(szs, gaps, "o-", color=PURP, lw=2)
    ax.set_xlabel("subgraph size"); ax.set_ylabel("greedy gap vs ILP (%)")
    ax.axhline(ilp["overall"]["mean_gap_pct"], color=GOLD, ls="--",
               label=f"mean {ilp['overall']['mean_gap_pct']}%")
    ax.legend(); style(ax, "D  ILP optimality gap (PuLP/CBC)")

    # E: centrality
    ax = fig.add_subplot(gs[1, 1])
    ax.bar(["circuit", "non-circuit"],
           [cent["circuit_mean"], cent["noncircuit_mean"]],
           color=[GOLD, BLUE], alpha=0.85)
    ax.set_ylabel("mean betweenness")
    style(ax, f"E  Centrality (lower in circuit, p={cent['p_value']:.3f})")

    # F: confidence tiers
    ax = fig.add_subplot(gs[1, 2])
    ax.plot([t["pct"] for t in tiers], [t["N"] for t in tiers], "o-",
            color=GREEN, lw=2)
    ax.set_xlabel("top-k% by NBLAST confidence"); ax.set_ylabel("MCIS N")
    style(ax, "F  Stability across confidence tiers")

    try:
        ac, an = annotation_quality()
        fig.text(0.5, 0.005,
                 f"Annotation (manual_cluster): circuit {ac:.1f}% vs non-circuit {an:.1f}%",
                 ha="center", fontsize=10)
    except Exception as e:
        print(f"(annotation panel skipped: {e})")

    fig.suptitle("MCIS Robustness · Statistical Validation · Optimality",
                 fontsize=15, fontweight="bold", y=0.99)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
